#[derive(Debug)]
pub enum RootError {
    ZeroDerivative(f64),
    MaxIterations,
}

/// Root finding using Newton-Raphson.  The 'f' and 'fder' are the function and 
/// its first derivative while 'start' indicates the initial guess.
pub fn newton_raphson<F1, F2>(
    f: &F1,
    fder: &F2,
    start: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    let epsilon = 1e-9;

    let mut x_pre = start;
    let mut x_cur = nr_iteration(f, fder, x_pre)?;
    let mut it = 1;

    loop {
        if (x_cur - x_pre).abs() < epsilon {
            return Ok(x_cur);
        }

        if it > max_iter {
            return Err(RootError::MaxIterations);
        }

        x_pre = x_cur;
        x_cur = nr_iteration(f, fder, x_pre)?;
        it += 1;
    }
}

/// Evaluate a single iteration for Newton's method.  Returns an error if the
/// derivative evaluates to zero.
fn nr_iteration<F1, F2>(f: &F1, fder: &F2, x: f64) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    let denom = fder(x);
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
        let fder = |_| 0.0;
        match newton_raphson(&f, &fder, 5.8, 100).expect_err("zero derivative not ok") {
            RootError::ZeroDerivative(_) => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    #[test]
    fn it_works() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let fder = |x| 2.0 * x - 9.0;

        let root = newton_raphson(&f, &fder, 5.8, 100).expect("found root");
        assert!((root - 5.0).abs() < 1e-9, "wanted root x=5");

        let root = newton_raphson(&f, &fder, 3.8, 100).expect("found root");
        assert!((root - 4.0).abs() < 1e-9, "wanted root x=4");
    }
}
