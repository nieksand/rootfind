/// Bounds represents the closed interval [a,b].
#[derive(Clone, Debug, PartialEq)]
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

    pub fn middle(&self) -> f64 {
        self.a + (self.b - self.a) / 2.0
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
    Singularity { x: f64 },
    IterationLimit { last_x: f64 },
}

/// Root finding via Bisection Method.
pub fn bisection<F>(f: &F, bounds: &Bounds, max_iter: usize) -> Result<f64, RootError>
where
    F: Fn(f64) -> f64,
{
    let mut window: Bounds = (*bounds).clone();
    let mut f_a = f(window.a);

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
    Err(RootError::IterationLimit { last_x: window.a })
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
    assert!(start.is_finite());

    let epsilon = 1e-9;

    let mut x_pre = start;
    let mut x_cur = nr_iteration(f, df, x_pre)?;
    let mut it = 1;

    loop {
        // convergence criteria
        if (x_cur - x_pre).abs() < epsilon {
            // sanity check
            if f(x_cur).abs() > epsilon {
                return Err(RootError::Singularity { x: x_cur });
            }
            return Ok(x_cur);
        }

        // stopping criteria
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
    #[should_panic]
    fn test_newton_nonfinite_start() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;
        let _ = newton_raphson(&f, &df, std::f64::NAN, 100);
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
    #[should_panic]
    fn test_bisection_no_straddle() {
        let f = |x| x * x;
        let _ = bisection(&f, &Bounds::new(-10.0, -5.0), 100);
    }

    #[test]
    fn test_bisection_parabola() {
        let f = |x| (x - 5.0) * (x - 4.0);

        let root = bisection(&f, &Bounds::new(4.5, 100.0), 100).expect("found root");
        assert!((root - 5.0).abs() < 1e-9, "wanted root x=5");

        let root = bisection(&f, &Bounds::new(-100000.0, 4.01), 100).expect("found root");
        assert!((root - 4.0).abs() < 1e-9, "wanted root x=4");
    }

    #[test]
    fn test_bisection_wikipedia() {
        // first example from wikipedia
        let f = |x| x * x - 612.0;
        let root = bisection(&f, &Bounds::new(10.0, 30.0), 100).expect("found root");
        assert!((root - 24.7386337537).abs() < 1e-9);

        // second example from wikipedia
        let f = |x: f64| x.cos() - x * x * x;
        let root = bisection(&f, &Bounds::new(0.0, 1.0), 100).expect("found root");
        assert!((root - 0.865474033102).abs() < 1e-9);
    }
}
