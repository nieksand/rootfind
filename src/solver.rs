//! Root finding algorithms.
//!
//! Fucntions typically have to be wrapped before use.  See the `wrap` module
//! for how to do this.
//!
//! Custom convergence criteria can be supplied.  Canned ones exist in the
//! `convergence` module.
//!
//! # Examples
//! Using Newton-Raphson:
//!
//! ```
//! use rootfind::solver::newton_raphson;
//! use rootfind::convergence::DeltaX;
//! use rootfind::wrap::RealFnAndFirst;
//!
//! // function and its derivative
//! let in_f = |x: f64| -x*x + 2.0*x + 1.0;
//! let in_df = |x: f64| -2.0*x + 2.0;
//! let f = RealFnAndFirst::new(&in_f, &in_df);
//!
//! // convergence criteria
//! let conv = DeltaX::new(1e-9);
//!
//! // invoke Newton-Raphson
//! let max_iterations = 20;
//! let root = newton_raphson(&f, 3.0, &conv, max_iterations).expect("root");
//!
//! // root at x=1+sqrt(2)
//! assert!((root-2.41421356237).abs() < 1e-9);
//! ```
//!
//! Using Bisection Method:
//!
//! ```
//! use rootfind::bracket::Bounds;
//! use rootfind::solver::bisection;
//!
//! // function... no derivatives needed!
//! let in_f = |x: f64| -x*x + 2.0*x + 1.0;
//!
//! // invoke bisection
//! let max_iterations = 20;
//! let root = bisection(&in_f, &Bounds::new(2.0, 3.0), 100).expect("root");
//!
//! // root at x=1+sqrt(2)
//! assert!((root-2.41421356237).abs() < 1e-9);
//! ```

use std::f64;
use bracket::{is_sign_change, Bounds};
use wrap::{RealD2fEval, RealDfEval, RealFnEval};
use convergence::IsConverged;

/// Root finding error conditions.
///
/// To help with diagnostics, these errors typically return the last relevant
/// `x` position.
#[derive(Debug)]
pub enum RootError {
    /// Derivative went to zero for method that depends on it to determine next
    /// step.
    ZeroDerivative { x: f64 },

    /// The solver computed a NaN for its next step x-value.
    IteratedToNaN { x: f64 },

    /// Iteration limit was reached.
    IterationLimit { last_x: f64 },
}

/// Driver for iterative root finders.
///
/// Allows for arbitrary iteration functions and converge criteria.  The user
/// function 'f' is kept compatible with the iteration routine using trait
/// bounds defined in 'wrap' module.
fn iterative_root_find<F, I, C>(
    f: &F,
    iterate: &I,
    start: f64,
    finish: &C,
    max_iter: usize,
) -> Result<f64, RootError>
where
    I: Fn(&F, f64) -> Result<(f64, f64), RootError>,
    C: IsConverged,
{
    assert!(start.is_finite());

    let mut x_pre = start;
    let mut x_cur = start;

    // stay inside maximum iteration count
    for _ in 0..max_iter {
        // invoke iteration method
        let t = iterate(f, x_pre)?;
        x_cur = t.0;
        let f_cur = t.1;

        // check convergence
        if finish.is_converged(x_pre, x_cur, f_cur) {
            return Ok(x_cur);
        }

        x_pre = x_cur;
    }
    return Err(RootError::IterationLimit { last_x: x_cur });
}

/// Root finding using Newton-Raphson.
///
/// The `start` indicates the initial guess.  For guesses sufficiently close to
/// the root this algorithm has quadratic convergence.
///
/// This algorithm requires the first derivative of f(x).
///
/// * If the second derivative is also available, consider Halley's method.
/// * If analytically computed derivatives are not available, consider Brent-Decker.
pub fn newton_raphson<F, C>(
    f: &F,
    start: f64,
    finish: &C,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F: RealFnEval + RealDfEval,
    C: IsConverged,
{
    iterative_root_find(f, &nr_iteration, start, finish, max_iter)
}

/// Evaluate a single iteration for Newton's method.  Returns an error if the
/// derivative evaluates to zero.  Returns (x_new, f(x_new)) otherwise.
fn nr_iteration<F>(f: &F, x: f64) -> Result<(f64, f64), RootError>
where
    F: RealFnEval + RealDfEval,
{
    let denom = f.eval_df(x);
    if denom == 0.0 {
        return Err(RootError::ZeroDerivative { x });
    }
    let f_x = f.eval_f(x);
    let x_new = x - f_x / denom;
    if !x_new.is_finite() {
        return Err(RootError::IteratedToNaN { x });
    }
    Ok((x_new, f_x))
}

/// Root finding using Halley's method.
///
/// The `start` indicates the initial guess.  For guesses sufficiently close to
/// the root this algorithm has cubic convergence.
///
/// This algorithm requires both the first and second derivatives of f(x).
///
/// * If only the first derivative is available, consider Newton-Raphson.
/// * If analytically computed derivatives are not available, consider Brent-Decker.
///
/// A good overview of the derivation, history, and geometric interpretation of
/// Halley's method is in:
///
/// *Scavo, T. R.; Thoo, J. B. (1995). "On the geometry of Halley's method".
/// American Mathematical Monthly. 102 (5): 417–426.*
///
pub fn halley_method<F, C>(f: &F, start: f64, finish: &C, max_iter: usize) -> Result<f64, RootError>
where
    F: RealFnEval + RealDfEval + RealD2fEval,
    C: IsConverged,
{
    iterative_root_find(f, &halley_iteration, start, finish, max_iter)
}

/// Evaluate a single iteration for Halley's method.  Returns (x_new, f(x_new))
/// on success.
fn halley_iteration<F>(f: &F, x: f64) -> Result<(f64, f64), RootError>
where
    F: RealFnEval + RealDfEval + RealD2fEval,
{
    let f_x = f.eval_f(x);
    let df_x = f.eval_df(x);
    let d2f_x = f.eval_d2f(x);

    if df_x == 0.0 {
        return Err(RootError::ZeroDerivative { x });
    }

    let x_new = x - (2.0 * f_x * df_x) / (2.0 * df_x * df_x - f_x * d2f_x);
    if !x_new.is_finite() {
        return Err(RootError::IteratedToNaN { x });
    }
    Ok((x_new, f_x))
}

/// Root finding via Bisection Method.
///
/// It always converges given a valid starting bracket, but the speed of
/// convergence is linear.
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

#[cfg(test)]
mod tests {
    use super::*;
    use convergence::DeltaX;
    use wrap::{RealFnAndFirst, RealFnAndFirstSecond};

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
                f: |x| x * x * x - x - 2.0,
                df: |x| 3.0 * x * x - 1.0,
                d2f: |x| 6.0 * x,
                roots: vec![1.52137970680457],
                guesses: vec![1.0],
                brackets: vec![Bounds::new(1.0, 2.0)],
            },
            RootTest {
                name: "Isaac Newton's Secant Example",
                f: |x| x * x * x + 10.0 * x * x - 7.0 * x - 44.0,
                df: |x| 3.0 * x * x + 20.0 * x - 7.0,
                d2f: |x| 6.0 * x + 20.0,
                roots: vec![2.20681731724844],
                guesses: vec![2.2],
                brackets: vec![Bounds::new(2.0, 2.3)],
            },
            RootTest {
                name: "Isaac Newton's NR Example",
                f: |x| x * x * x - 2.0 * x - 5.0,
                df: |x| 3.0 * x * x - 2.0,
                d2f: |x| 6.0 * x,
                roots: vec![2.0945514815423265],
                guesses: vec![2.0],
                brackets: vec![Bounds::new(2.0, 3.0)],
            },
            RootTest {
                name: "Thomas Simpson NR Example",
                f: |x| {
                    (1. - x).sqrt() + (1. - 2. * x * x).sqrt() + (1. - 3. * x * x * x).sqrt() - 2.
                },
                df: |x| {
                    -2. * x * (1. - 2. * x * x).sqrt().recip()
                        - 9. * x * x * (1. - 3. * x * x * x).sqrt().recip() / 2.
                        - 1. * (1. - x).sqrt().recip() / 2.
                },
                d2f: |x| {
                    -9. * x * (1. - 3. * x * x * x).sqrt().recip()
                        - 4. * x * x * (1. - 2. * x * x).powf(1.5)
                        - 2. * (1. - 2. * x * x).sqrt().recip()
                        - 81. * x * x * x * x * (1. - 3. * x * x * x).powf(1.5) / 4.
                        - 1. * (1. - x).powf(1.5) / 4.
                },
                roots: vec![0.55158615249704711724768527],
                guesses: vec![0.5],
                brackets: vec![Bounds::new(0.0, 3.0)],
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
        let conv = DeltaX::new(1e-9);
        for t in make_root_tests() {
            for i in 0..t.roots.len() {
                let f = RealFnAndFirst::new(&t.f, &t.df);
                let root = newton_raphson(&f, t.guesses[i], &conv, 100).expect("found root");
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
        let in_f = |x| (x - 5.0) * (x - 4.0);
        let in_df = |x| 2.0 * x - 9.0;
        let f = RealFnAndFirst::new(&in_f, &in_df);

        let conv = DeltaX::new(1e-9);
        let _ = newton_raphson(&f, f64::NAN, &conv, 100);
    }

    #[test]
    fn test_newton_zero_derivative() {
        let in_f = |_| 2.0;
        let in_df = |_| 0.0;
        let f = RealFnAndFirst::new(&in_f, &in_df);

        let conv = DeltaX::new(1e-9);
        match newton_raphson(&f, 5.8, &conv, 100).expect_err("zero derivative not ok") {
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
        let conv = DeltaX::new(1e-9);
        for t in make_root_tests() {
            for i in 0..t.roots.len() {
                let f = RealFnAndFirstSecond::new(&t.f, &t.df, &t.d2f);
                let root = halley_method(&f, t.guesses[i], &conv, 100).expect("found root");
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
        let in_f = |x: f64| x.sin();
        let in_df = |x: f64| x.cos();
        let in_d2f = |x: f64| -x.sin();
        let f = RealFnAndFirstSecond::new(&in_f, &in_df, &in_d2f);

        let conv = DeltaX::new(1e-9);
        let _ = halley_method(&f, f64::NAN, &conv, 100);
    }

    #[test]
    fn test_pathology_microstep() {
        // f(x) = 0.01*e^(1/x)-1
        let in_f = |x: f64| 0.001 * (1.0 / x).exp() - 1.0;
        let in_df = |x: f64| -0.001 * (1.0 / x).exp() / (x * x);
        let f = RealFnAndFirst::new(&in_f, &in_df);

        let conv = DeltaX::new(1e-9);
        let root = newton_raphson(&f, 0.00142, &conv, 100).expect("root");

        // we "converged" but are far from actual root
        assert!(root.abs() > 0.001);
        assert!((root - 0.144765).abs() > 0.14);
    }

    #[test]
    fn test_pathology_flatlining() {
        // f(x)=1/e^(x^100) - 0.5, roots approx -0.996342 and 0.996342
        let in_f = |x: f64| x.powi(-100).exp() - 0.5;

        // d/dx = -100e^(-x^100) * x^99
        let in_df = |x: f64| -100.0 * (-x.powi(100)).exp() * x.powi(99);

        let f = RealFnAndFirst::new(&in_f, &in_df);

        let conv = DeltaX::new(1e-9);
        let _ = newton_raphson(&f, 0.99999, &conv, 100).expect_err("no convergence");
    }
}
