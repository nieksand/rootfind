/// Bounds represents interval [a,b).
#[derive(Debug, PartialEq)]
pub struct Bounds {
    a: f64,
    b: f64,
}

impl Bounds {
    pub fn new(a: f64, b: f64) -> Bounds {
        assert!(a < b);
        assert!(a.is_finite() && b.is_finite());
        Bounds { a, b }
    }
}

/// Use a sliding window technique for find an interval where the signs change
/// on the extents.  If f is continuous, then the Intermediate Value Theorem
/// guarantees our interval contains at least one root.
pub fn sign_change_window<F>(f: &F, bounds: &Bounds, window_size: f64) -> Option<Bounds>
where
    F: Fn(f64) -> f64,
{
    let mut win = Bounds {
        a: bounds.a,
        b: bounds.a + window_size,
    };

    let mut f_a = f(win.a);
    while win.a < bounds.b {
        let f_b = f(win.b);

        // sign change at current window
        if f_a == 0.0 || f_a * f_b < 0.0 {
            return Some(win);
        }

        f_a = f_b;
        win.a = win.b;
        win.b += window_size;
    }
    None
}

#[derive(Debug)]
pub enum RootError {
    ZeroDerivative(f64),
    MaxIterations,
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
            return Err(RootError::MaxIterations);
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
        return Err(RootError::ZeroDerivative(x));
    }
    Ok(x - f(x) / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_change_window_hit() {
        // root at x=-9
        let f = |x| x + 9.0;
        let win = sign_change_window(&f, &Bounds::new(-100.0, 100.0), 10.0).expect("window found");
        assert_eq!(win, Bounds::new(-10.0, 0.0));

        // sign change on inclusive-side of window boundary
        let win = sign_change_window(&f, &Bounds::new(-29.0, -8.0), 10.0).expect("window found");
        assert_eq!(win, Bounds::new(-9.0, 1.0));
    }

    #[test]
    fn test_sign_change_window_miss() {
        // root at x=-9, but window doesn't include
        let f = |x| x + 9.0;
        let win = sign_change_window(&f, &Bounds::new(0.0, 100.0), 10.0);
        assert!(win.is_none());

        // sign change on exclusive-side of window boundary
        let win = sign_change_window(&f, &Bounds::new(-29.0, -9.0), 10.0);
        assert!(win.is_none());

        // no root
        let f = |_| 33.0;
        let win = sign_change_window(&f, &Bounds::new(-100.0, 100.0), 1.0);
        assert!(win.is_none());
    }

    #[test]
    fn test_zero_derivative() {
        let f = |_| 2.0;
        let df = |_| 0.0;
        match newton_raphson(&f, &df, 5.8, 100).expect_err("zero derivative not ok") {
            RootError::ZeroDerivative(_) => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    #[test]
    fn test_parabola() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;

        let root = newton_raphson(&f, &df, 5.8, 100).expect("found root");
        assert!((root - 5.0).abs() < 1e-9, "wanted root x=5");

        let root = newton_raphson(&f, &df, 3.8, 100).expect("found root");
        assert!((root - 4.0).abs() < 1e-9, "wanted root x=4");
    }

    #[test]
    fn test_wikipedia() {
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
}
