use super::RootError;
use bracket::{is_sign_change, Bounds};
use convergence::IsConverged;
use wrap::RealFnEval;

/// Driver for iterative root finders.
///
/// Allows for arbitrary iteration functions and converge criteria.  The user
/// function 'f' is kept compatible with the iteration routine using trait
/// bounds defined in 'wrap' module.
pub fn iterative_root_find<F, I, C>(
    f: &F,
    iterate: &I,
    start: f64,
    finish: &C,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F: RealFnEval,
    I: Fn(&F, f64, f64) -> Result<f64, RootError>,
    C: IsConverged,
{
    assert!(start.is_finite());

    let mut x_pre = start;
    let mut x_cur = start;

    let mut f_pre = f.eval_f(x_pre);
    let mut f_cur;

    // stay inside maximum iteration count
    for _ in 0..max_iter {
        // invoke iteration method
        x_cur = iterate(f, x_pre, f_pre)?;
        f_cur = f.eval_f(x_cur);

        // check convergence
        if finish.is_converged(x_pre, x_cur, f_cur) {
            return Ok(x_cur);
        }

        x_pre = x_cur;
        f_pre = f_cur;
    }
    return Err(RootError::IterationLimit { last_x: x_cur });
}

/// Safe solver hybidizes iterative method to ensure convergence.
///
/// Do not use this!
///
/// The code is both messy and untrustworthy.  I had to tighten to convergence
/// criteria to get expected results, for reasons I don't yet understand.
///
pub fn safe_iterative_root_find<F, I, C>(
    f: &F,
    iterate: &I,
    bounds: &Bounds,
    finish: &C,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F: RealFnEval,
    I: Fn(&F, f64, f64) -> Result<f64, RootError>,
    C: IsConverged,
{
    // solver always kept in window
    let mut window = bounds.clone();
    let mut f_a = f.eval_f(window.a);
    assert!(is_sign_change(f_a, f.eval_f(window.b)));

    // bisect once to start in middle of bracket
    let mut must_bisect = true;
    let mut x_pre = window.a;
    let mut f_pre = f_a;

    // previous window sizes
    let mut sz_pre = window.size();
    let mut sz_pre2 = window.size();

    println!("\n== SOLVER");
    for _ in 0..max_iter {
        let x_cur;
        let f_cur;

        if must_bisect {
            println!("--> bisection");
            x_cur = window.middle();
            f_cur = f.eval_f(x_cur);

            if is_sign_change(f_a, f_cur) {
                window.b = x_cur;
            } else {
                window.a = x_cur;
                f_a = f_cur;
            }
            must_bisect = false;
        } else {
            match iterate(f, x_pre, f_pre) {
                Ok(x_new) if window.contains(x_new) => {
                    println!("--> iteration");
                    x_cur = x_new;
                    f_cur = f.eval_f(x_cur);

                    if is_sign_change(f_a, f_cur) {
                        window.b = x_cur;
                    } else {
                        window.a = x_cur;
                        f_a = f_cur;
                    }
                }
                _ => {
                    println!("--> failed iteration");
                    // out of bounds or other error (e.g. zero derivative)
                    must_bisect = true;
                    continue;
                }
            }
        }

        println!("~~> [{}, {}]", window.a, window.b);
        // check convergence
        //if window.size() < 1e-9 {
        //    return Ok(window.middle());
        //}
        if finish.is_converged(x_pre, x_cur, f_cur) {
            return Ok(x_cur);
        }

        // insufficient progress
        let sz_cur = window.size();
        if sz_pre2 * 0.5 < sz_cur {
            println!("--> SHRINK YOU APE!");
            must_bisect = true;
        }

        sz_pre2 = sz_pre;
        sz_pre = sz_cur;
        x_pre = x_cur;
        f_pre = f_cur;
    }

    Err(RootError::IterationLimit {
        last_x: window.middle(),
    })
}
