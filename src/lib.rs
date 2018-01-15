/// Bounds represents the closed interval [a,b].
#[derive(Debug, PartialEq)]
pub struct Bounds {
    a: f64,
    b: f64,
}

impl Bounds {
    pub fn new(a: f64, b: f64) -> Bounds {
        assert!(a <= b);
        assert!(a.is_finite() && b.is_finite());
        Bounds { a, b }
    }
}

/// Whether signs of values differ, properly handling integer underflow.
fn is_sign_change(lhs: f64, rhs: f64) -> bool {
    lhs.signum() != rhs.signum()
}

/// Scans interval [a,b] and emits the first bracket containing a sign change.
/// For a continuous function the Intermediate Value Theorem guarantees that the
/// bracket contains at least one root.  Without a continuity guarantee, it
/// might be a singularity instead.
pub fn first_bracket<F>(f: &F, bounds: &Bounds, window_size: f64) -> Option<Bounds>
where
    F: Fn(f64) -> f64,
{
    assert!(window_size > 0.0);

    let mut win = Bounds {
        a: bounds.a,
        b: (bounds.a + window_size).min(bounds.b),
    };

    let mut f_a = f(win.a);
    while win.a < bounds.b {
        let f_b = f(win.b);

        // found root or singularity
        if is_sign_change(f_a, f_b) {
            return Some(win);
        }

        f_a = f_b;
        win.a = win.b;
        win.b = (win.b + window_size).min(bounds.b);
    }
    None
}

#[derive(Debug)]
pub enum RootError {
    ZeroDerivative { x: f64 },
    IterationLimit { last_x: f64 },
}

/// Root finding using Newton-Raphson.  The 'f' and 'df' are the function and
/// its first derivative while 'start' indicates the initial guess.
pub fn newton_raphson<F1, F2>(
    f: &F1,
    df: &F2,
    start: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    let epsilon = 1e-9;

    let mut x_pre = start;
    let mut x_cur = nr_iteration(f, df, x_pre)?;
    let mut it = 1;

    loop {
        if (x_cur - x_pre).abs() < epsilon {
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
    Ok(x - f(x) / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bounds_new_valid() {
        let b = Bounds::new(-2.0, 2.0);
        assert_eq!(b.a, -2.0);
        assert_eq!(b.b, 2.0);

        let b = Bounds::new(2.0, 2.0);
        assert_eq!(b.a, 2.0);
        assert_eq!(b.b, 2.0);
    }

    #[test]
    #[should_panic]
    fn test_bounds_new_flipped_extents() {
        Bounds::new(2.0, -2.0);
    }

    #[test]
    #[should_panic]
    fn test_bounds_new_nan() {
        Bounds::new(std::f64::NAN, -2.0);
    }

    #[test]
    #[should_panic]
    fn test_bounds_new_infinite() {
        Bounds::new(std::f64::NEG_INFINITY, std::f64::INFINITY);
    }

    #[test]
    fn test_is_sign_change() {
        // easy peasy
        assert_eq!(is_sign_change(-1.0, -1.0), false);
        assert_eq!(is_sign_change(1.0, 1.0), false);
        assert_eq!(is_sign_change(-1.0, 1.0), true);

        // zero tests
        assert_eq!(is_sign_change(0.0, 0.0), false);
        assert_eq!(is_sign_change(0.0, 1.0), false);
        assert_eq!(is_sign_change(0.0, -1.0), true);

        // naughty signed zeroes
        assert_eq!(is_sign_change(-0.0, -1.0), false);
        assert_eq!(is_sign_change(-0.0, 0.0), true);
    }

    #[test]
    fn test_is_sign_change_underflow() {
        // floating point underflow breaks naive a*b<0 check
        assert_eq!(
            is_sign_change(1e-120, -2e-300),
            true,
            "sign change with float underflow"
        );
    }

    #[test]
    #[should_panic]
    fn test_first_bracket_negative_window() {
        let f = |x| x * x;
        first_bracket(&f, &Bounds::new(-20.0, 20.0), -1.0);
    }

    #[test]
    #[should_panic]
    fn test_first_bracket_zero_window() {
        let f = |x| x * x;
        first_bracket(&f, &Bounds::new(-20.0, 20.0), 0.0);
    }

    #[test]
    fn test_first_bracket_hit() {
        // root at x=-9
        let f = |x| x + 9.0;
        let win = first_bracket(&f, &Bounds::new(-100.0, 100.0), 10.0).expect("window found");
        assert_eq!(win, Bounds::new(-10.0, 0.0));

        // sign change on right window boundary
        let win = first_bracket(&f, &Bounds::new(-29.0, -8.0), 10.0).expect("window found");
        assert_eq!(win, Bounds::new(-19.0, -9.0));

        // sign change on left window boundary
        let win = first_bracket(&f, &Bounds::new(-19.0, -9.0), 10.0).expect("window found");
        assert_eq!(win, Bounds::new(-19.0, -9.0));
    }

    #[test]
    fn test_first_bracket_miss() {
        // root at x=-9, but window doesn't include
        let f = |x| x + 9.0;
        let win = first_bracket(&f, &Bounds::new(0.0, 100.0), 10.0);
        assert!(win.is_none());

        // no root
        let f = |_| 33.0;
        let win = first_bracket(&f, &Bounds::new(-100.0, 100.0), 1.0);
        assert!(win.is_none());
    }

    #[test]
    fn test_newton_zero_derivative() {
        let f = |_| 2.0;
        let df = |_| 0.0;
        match newton_raphson(&f, &df, 5.8, 100).expect_err("zero derivative not ok") {
            RootError::ZeroDerivative { .. } => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    #[test]
    fn test_newton_parabola() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;

        let root = newton_raphson(&f, &df, 5.8, 100).expect("found root");
        assert!((root - 5.0).abs() < 1e-9, "wanted root x=5");

        let root = newton_raphson(&f, &df, 3.8, 100).expect("found root");
        assert!((root - 4.0).abs() < 1e-9, "wanted root x=4");
    }

    #[test]
    fn test_newton_wikipedia() {
        // first example from wikipedia
        let f = |x| x * x - 612.0;
        let df = |x| 2.0 * x;
        let root = newton_raphson(&f, &df, 10.0, 100).expect("found root");
        assert!((root - 24.7386337537).abs() < 1e-9);

        // second example from wikipedia
        let f = |x: f64| x.cos() - x * x * x;
        let df = |x: f64| -x.sin() - 3.0 * x * x;
        let root = newton_raphson(&f, &df, 0.5, 100).expect("found root");
        assert!((root - 0.865474033102).abs() < 1e-9);
    }

    #[test]
    fn test_newton_singularity() {
        let f = |x| if x < 1.0 { x - 2.0 } else { x };
        let df = |x| if x != 1.0 { 1.0 } else { std::f64::NAN };
        let _sing = newton_raphson(&f, &df, 0.9, 100).expect("potato");
    }
}
