use std::f64;
use bracket::{is_sign_change, Bounds};

#[derive(Debug)]
pub enum RootError {
    ZeroDerivative { x: f64 },
    IteratedToNaN { x: f64 },
    ConvergedOnNonZero { x: f64 },
    IterationLimit { last_x: f64 },
}

/// Root finding via Bisection Method.  It converges globally but is relatively
/// slow.
pub fn bisection<F>(f: &F, bounds: &Bounds, max_iter: usize) -> Result<f64, RootError>
where
    F: Fn(f64) -> f64,
{
    let mut window: Bounds = (*bounds).clone();
    let mut f_a = f(window.a);

    // ensure we started with valid bracket
    assert!(is_sign_change(f_a, f(window.b)));

    for _ in 0..max_iter {
        let mid = window.middle();
        let f_mid = f(mid);

        if is_sign_change(f_a, f_mid) {
            window.b = mid;
        } else {
            window.a = mid;
            f_a = f_mid;
        }

        // convergence criteria
        if window.b - window.a < 1e-9 {
            return Ok(window.a);
        }
    }
    Err(RootError::IterationLimit {
        last_x: window.middle(),
    })
}

/// Root finding using Newton-Raphson.  For guesses sufficiently close to the
/// root this algorithm has quadratic convergence.
///
/// This algorithm requires the first derivative of f(x).  If the second
/// derivative is also available, consider Halley's method.  If no analytically
/// computed derivatives are available, consider Brent-Decker.
///
/// The 'f' and 'df' are the function and
/// its first derivative while 'start' indicates the initial guess.  The
/// 'accuracy' bounds how much x_old and x_new may differ and the max
/// |f(x_final|) before we declare convergence.
pub fn newton_raphson<F1, F2>(
    f: &F1,
    df: &F2,
    start: f64,
    accuracy: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    assert!(start.is_finite());
    assert!(accuracy > 0.0);

    let mut x_pre = start;
    let mut x_cur = nr_iteration(f, df, x_pre)?;
    let mut it = 1;

    loop {
        // convergence criteria
        if (x_cur - x_pre).abs() < accuracy {
            // possible if df is huge
            if f(x_cur) > accuracy {
                return Err(RootError::ConvergedOnNonZero { x: x_cur });
            }
            return Ok(x_cur);
        }

        if it > max_iter {
            return Err(RootError::IterationLimit { last_x: x_cur });
        }
        x_pre = x_cur;
        x_cur = nr_iteration(f, df, x_pre)?;
        it += 1;
    }
}

/// Evaluate a single iteration for Newton's method.  Returns an error if the
/// derivative evaluates to zero.
fn nr_iteration<F1, F2>(f: &F1, df: &F2, x: f64) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    let denom = df(x);
    if denom == 0.0 {
        return Err(RootError::ZeroDerivative { x });
    }
    let x_new = x - f(x) / denom;
    if !x_new.is_finite() {
        return Err(RootError::IteratedToNaN { x });
    }
    Ok(x_new)
}

/// Root finding using Halley's method.  For guesses sufficiently close to the
/// root this algorithm has cubic convergence.
///
/// This algorithm requires both the first and second derivatives of f(x).  If
/// only the first derivative is available, consider Newton-Raphson.  If no
/// analytically computed derivatives are available, consider Brent-Decker.
///
/// The 'f', 'df', and 'd2f' are the function and its first and second
/// derivatives.  The 'start' indicates the initial guess.  The 'accuracy'
/// bounds how much x_old and x_new may differ and the max |f(x_final|) before
/// we declare convergence.
///
/// A good overview of the derivation, history, and geometric interpretation of
/// Halley's method is in:
///
/// Scavo, T. R.; Thoo, J. B. (1995). "On the geometry of Halley's method".
/// American Mathematical Monthly. 102 (5): 417–426.
///
pub fn halley_method<F1, F2, F3>(
    f: &F1,
    df: &F2,
    d2f: &F3,
    start: f64,
    accuracy: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
    F3: Fn(f64) -> f64,
{
    assert!(start.is_finite());
    assert!(accuracy > 0.0);

    let mut x_pre = start;
    let mut x_cur = halley_iteration(f, df, d2f, x_pre)?;
    let mut it = 1;

    loop {
        // convergence criteria
        if (x_cur - x_pre).abs() < accuracy {
            // possible if df is huge and d2f near zero
            if f(x_cur) > accuracy {
                return Err(RootError::ConvergedOnNonZero { x: x_cur });
            }
            return Ok(x_cur);
        }

        if it > max_iter {
            return Err(RootError::IterationLimit { last_x: x_cur });
        }
        x_pre = x_cur;
        x_cur = halley_iteration(f, df, d2f, x_pre)?;
        it += 1;
    }
}

/// Evaluate a single iteration for Halley's method.
fn halley_iteration<F1, F2, F3>(f: &F1, df: &F2, d2f: &F3, x: f64) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
    F3: Fn(f64) -> f64,
{
    let f_x = f(x);
    let df_x = df(x);
    let d2f_x = d2f(x);

    if df_x == 0.0 {
        return Err(RootError::ZeroDerivative { x });
    }

    let x_new = x - (2.0 * f_x * df_x) / (2.0 * df_x * df_x - f_x * d2f_x);
    if !x_new.is_finite() {
        return Err(RootError::IteratedToNaN { x });
    }
    Ok(x_new)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct RootTest {
        name: &'static str,
        f: fn(f64) -> f64,
        df: fn(f64) -> f64,
        d2f: fn(f64) -> f64,
        roots: Vec<f64>,
        guesses: Vec<f64>,
        brackets: Vec<Bounds>,
    }

    fn make_root_tests() -> Vec<RootTest> {
        vec![
            RootTest {
                name: "Factored Parabola",
                f: |x| (x - 5.0) * (x - 4.0),
                df: |x| 2.0 * x - 9.0,
                d2f: |_| 2.0,
                roots: vec![5.0, 4.0],
                guesses: vec![5.8, 3.8],
                brackets: vec![Bounds::new(4.5, 100.0), Bounds::new(-100000.0, 4.01)],
            },
            RootTest {
                name: "Wikipedia NR Parabola",
                f: |x| x * x - 612.0,
                df: |x| 2.0 * x,
                d2f: |_| 2.0,
                roots: vec![-24.7386337537, 24.7386337537],
                guesses: vec![-10.0, 10.0],
                brackets: vec![Bounds::new(-30.0, 10.0), Bounds::new(10.0, 30.0)],
            },
            RootTest {
                name: "Wikipedia NR Trigonometry",
                f: |x| x.cos() - x * x * x,
                df: |x| -x.sin() - 3. * x * x,
                d2f: |x| -x.cos() - 6. * x,
                roots: vec![0.865474033102],
                guesses: vec![0.5],
                brackets: vec![Bounds::new(0.0, 1.0)],
            },
            RootTest {
                name: "Wikipedia Bisection Cubic",
                f: |x| x*x*x - x - 2.0,
                df: |x| 3.0*x*x - 1.0,
                d2f: |x| 6.0*x,
                roots: vec![1.52137970680457],
                guesses: vec![1.0],
                brackets: vec![Bounds::new(1.0, 2.0),],
            },
        ]
    }

    #[test]
    fn test_bisection_root_finding() {
        for t in make_root_tests() {
            for i in 0..t.roots.len() {
                let root = bisection(&t.f, &t.brackets[i], 100).expect("found root");
                assert!(
                    (root - t.roots[i]).abs() < 1e-9,
                    format!("{} root wanted={}, got={}", t.name, t.roots[i], root)
                );
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_bisection_no_straddle() {
        let f = |x| x * x;
        let _ = bisection(&f, &Bounds::new(-10.0, -5.0), 100);
    }

    #[test]
    fn test_bisection_centered_root() {
        let f = |x| x;
        let root = bisection(&f, &Bounds::new(-1000000.0, 1000000.0), 100).expect("found root");
        assert!(root.abs() < 1e-9, "wanted root x=0");
    }

    #[test]
    fn test_newton_root_finding() {
        for t in make_root_tests() {
            for i in 0..t.roots.len() {
                let root =
                    newton_raphson(&t.f, &t.df, t.guesses[i], 1e-9, 100).expect("found root");
                assert!(
                    (root - t.roots[i]).abs() < 1e-9,
                    format!("{} root wanted={}, got={}", t.name, t.roots[i], root)
                );
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_newton_nonfinite_start() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;
        let _ = newton_raphson(&f, &df, f64::NAN, 1e-9, 100);
    }

    #[test]
    #[should_panic]
    fn test_newton_accuracy_negative() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;
        let _ = newton_raphson(&f, &df, 42.0, -1e-9, 100);
    }

    #[test]
    #[should_panic]
    fn test_newton_accuracy_zero() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;
        let _ = newton_raphson(&f, &df, 42.0, 0.0, 100);
    }

    #[test]
    fn test_newton_zero_derivative() {
        let f = |_| 2.0;
        let df = |_| 0.0;
        match newton_raphson(&f, &df, 5.8, 1e-9, 100).expect_err("zero derivative not ok") {
            RootError::ZeroDerivative { .. } => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    #[test]
    fn test_halley_root_finding() {
        for t in make_root_tests() {
            for i in 0..t.roots.len() {
                let root = halley_method(&t.f, &t.df, &t.d2f, t.guesses[i], 1e-9, 100)
                    .expect("found root");
                assert!(
                    (root - t.roots[i]).abs() < 1e-9,
                    format!("{} root wanted={}, got={}", t.name, t.roots[i], root)
                );
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_halley_nonfinite_start() {
        let f = |x: f64| x.sin();
        let df = |x: f64| x.cos();
        let d2f = |x: f64| -x.sin();
        let _ = halley_method(&f, &df, &d2f, f64::NAN, 1e-9, 100);
    }

    #[test]
    #[should_panic]
    fn test_halley_accuracy_negative() {
        let f = |x: f64| x.sin();
        let df = |x: f64| x.cos();
        let d2f = |x: f64| -x.sin();
        let _ = halley_method(&f, &df, &d2f, 42.0, -1e-9, 100);
    }

    #[test]
    #[should_panic]
    fn test_halley_accuracy_zero() {
        let f = |x: f64| x.sin();
        let df = |x: f64| x.cos();
        let d2f = |x: f64| -x.sin();
        let _ = halley_method(&f, &df, &d2f, 42.0, 0.0, 100);
    }

    #[test]
    fn test_pathology_microstep() {
        // f(x) = 0.01*e^(1/x)-1
        let f = |x: f64| 0.001 * (1.0 / x).exp() - 1.0;
        let df = |x: f64| -0.001 * (1.0 / x).exp() / (x * x);
        match newton_raphson(&f, &df, 0.00142, 1e-9, 100).expect_err("microstep fail") {
            RootError::ConvergedOnNonZero { .. } => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    #[test]
    fn test_pathology_flatlining() {
        // f(x)=1/e^(x^100) - 0.5, roots approx -0.996342 and 0.996342
        let f = |x: f64| x.powi(-100).exp() - 0.5;

        // d/dx = -100e^(-x^100) * x^99
        let df = |x: f64| -100.0 * (-x.powi(100)).exp() * x.powi(99);

        let _ = newton_raphson(&f, &df, 0.99999, 1e-9, 100).expect_err("no convergence");
    }
}
