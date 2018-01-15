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
