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
        self.a + (self.b - self.a) * 0.5
    }
}

/// BracketGenerator iterates over all root-holding brackets.  Internally it is
/// making repeated calls to first_bracket until the entire bounds are explored.
pub struct BracketGenerator<'a, F: 'a> {
    f: &'a F,
    remaining: Option<Bounds>,
    window_size: f64,
}

impl<'a, F> BracketGenerator<'a, F>
where
    F: Fn(f64) -> f64,
{
    pub fn new(f: &F, bounds: Bounds, window_size: f64) -> BracketGenerator<F> {
        BracketGenerator {
            f,
            remaining: Some(bounds),
            window_size,
        }
    }
}

impl<'a, F> Iterator for BracketGenerator<'a, F>
where
    F: Fn(f64) -> f64,
{
    type Item = Bounds;

    fn next(&mut self) -> Option<Bounds> {
        let mut search_bounds = self.remaining.clone()?;
        let result = first_bracket(&self.f, &search_bounds, self.window_size);

        match result {
            None => {
                self.remaining = None;
            }
            Some(ref found_bracket) => {
                search_bounds.a = found_bracket.b;
                self.remaining = Some(search_bounds);
            }
        }
        result
    }
}

/// Whether signs differ, properly handling integer underflow.
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
    IteratedToNaN { x: f64 },
    ConvergedOnNonZero { x: f64 },
    IterationLimit { last_x: f64 },
}

/// Root finding via Bisection Method.
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
            // maybe first derivative is huge
            if f(x_cur) > epsilon {
                return Err(RootError::ConvergedOnNonZero { x: x_cur });
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
    let x_new = x - f(x) / denom;
    if !x_new.is_finite() {
        return Err(RootError::IteratedToNaN { x });
    }
    Ok(x_new)
}

/// Root finding using Halley's method.  The 'f', 'df', and 'd2f' are the
/// function and it's first and second derivatives.  The 'start' indicates the
/// initial guess.
pub fn halley_method<F1, F2, F3>(
    f: &F1,
    df: &F2,
    d2f: &F3,
    start: f64,
    max_iter: usize,
) -> Result<f64, RootError>
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
    F3: Fn(f64) -> f64,
{
    assert!(start.is_finite());

    let epsilon = 1e-9;

    let mut x_pre = start;
    let mut x_cur = halley_iteration(f, df, d2f, x_pre)?;
    let mut it = 1;

    loop {
        // convergence criteria
        if (x_cur - x_pre).abs() < epsilon {
            // maybe first derivative is huge and second derivitive near zero
            if f(x_cur) > epsilon {
                return Err(RootError::ConvergedOnNonZero { x: x_cur });
            }
            return Ok(x_cur);
        }

        // stopping criteria
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

    Ok(x - (2.0 * f_x * df_x) / (2.0 * df_x * df_x - f_x * d2f_x))
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
    fn test_bracket_generator_hits() {
        let f = |x: f64| x.sin();
        let pi = std::f64::consts::PI;
        let b = Bounds::new(-0.1, 4.0 * pi + 0.1);

        let results: Vec<Bounds> = BracketGenerator::new(&f, b, 1.0).collect();

        // should bracket 0, pi, 2pi, 3pi, 4pi
        assert_eq!(results.len(), 5);

        assert_eq!(Bounds::new(-0.1, 0.9), results[0]); // 0
        assert_eq!(Bounds::new(2.9, 3.9), results[1]); // pi
        assert_eq!(Bounds::new(5.9, 6.9), results[2]); // 2pi
        assert_eq!(Bounds::new(8.9, 9.9), results[3]); // 3pi
        assert_eq!(Bounds::new(11.9, 4.0 * pi + 0.1), results[4]); // 4pi
    }

    #[test]
    fn test_bracket_generator_empty() {
        let f = |x: f64| x.sin();
        let b = Bounds::new(0.1, 0.5);

        let mut gen = BracketGenerator::new(&f, b, 0.1);
        assert!(gen.next().is_none());
    }

    #[test]
    #[should_panic]
    fn test_bracket_generator_window_negative() {
        let f = |x: f64| x.sin();
        let b = Bounds::new(-10.0, 10.0);
        let _brackets: Vec<Bounds> = BracketGenerator::new(&f, b, -0.1).collect();
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

    struct RootTest<F1, F2, F3> {
        name: String,
        f: F1,
        df: F2,
        d2f: F3,
        roots: Vec<f64>,
        guesses: Vec<f64>,
        brackets: Vec<Bounds>,
    }

    type Rfn = Box<Fn(f64) -> f64>;

    /// Common set of tests for root-finding algorithms.  In order to satisfy
    /// various method requirements, it includes 1st and 2nd derivatives as well
    /// as (expected convergent) initial guesses and bounds.
    fn make_root_tests() -> Vec<RootTest<Rfn, Rfn, Rfn>> {
        vec![
            RootTest {
                name: "Factored Parabola".to_owned(),
                f: Box::new(|x| (x - 5.0) * (x - 4.0)),
                df: Box::new(|x| 2.0 * x - 9.0),
                d2f: Box::new(|_| 2.0),
                roots: vec![5.0, 4.0],
                guesses: vec![5.8, 3.8],
                brackets: vec![Bounds::new(4.5, 5.5), Bounds::new(3.5, 4.5)],
            },
            // first example from Wikipedia "Newton's Method" page
            RootTest {
                name: "Wikipedia Parabola".to_owned(),
                f: Box::new(|x| x * x - 612.0),
                df: Box::new(|x| 2.0 * x),
                d2f: Box::new(|_| 2.0),
                roots: vec![-24.7386337537, 24.7386337537],
                guesses: vec![-10.0, 10.0],
                brackets: vec![Bounds::new(-100.0, 0.0), Bounds::new(0.0, 500.0)],
            },
            // second example from Wikipedia "Newton's Method" page
            RootTest {
                name: "Wikipedia Trigonometry".to_owned(),
                f: Box::new(|x| x.cos() - x * x * x),
                df: Box::new(|x| -x.sin() - 3.0 * x * x),
                d2f: Box::new(|x| -x.cos() - 6.0 * x),
                roots: vec![0.865474033102],
                guesses: vec![0.5],
                brackets: vec![Bounds::new(0.0, 10.0)],
            },
        ]
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
    fn test_halley_parabola() {
        let f = |x| (x - 5.0) * (x - 4.0);
        let df = |x| 2.0 * x - 9.0;
        let d2f = |_| 2.0;

        let root = halley_method(&f, &df, &d2f, 5.8, 100).expect("found root");
        assert!((root - 5.0).abs() < 1e-9, "wanted root x=5");

        let root = halley_method(&f, &df, &d2f, 3.8, 100).expect("found root");
        assert!((root - 4.0).abs() < 1e-9, "wanted root x=4");
    }

    #[test]
    fn test_halley_wikipedia() {
        // first example from wikipedia
        let f = |x| x * x - 612.0;
        let df = |x| 2.0 * x;
        let d2f = |_| 2.0;
        let root = halley_method(&f, &df, &d2f, 10.0, 100).expect("found root");
        assert!((root - 24.7386337537).abs() < 1e-9);

        // second example from wikipedia
        let f = |x: f64| x.cos() - x * x * x;
        let df = |x: f64| -x.sin() - 3.0 * x * x;
        let d2f = |x: f64| -x.cos() - 6.0 * x;
        let root = halley_method(&f, &df, &d2f, 0.5, 100).expect("found root");
        assert!((root - 0.865474033102).abs() < 1e-9);
    }

    #[test]
    fn test_pathology_microstep() {
        // f(x) = 0.01*e^(1/x)-1
        let f = |x: f64| 0.001 * (1.0 / x).exp() - 1.0;
        let df = |x: f64| -0.001 * (1.0 / x).exp() / (x * x);
        match newton_raphson(&f, &df, 0.00142, 100).expect_err("microstep fail") {
            RootError::ConvergedOnNonZero { .. } => {
                return;
            }
            _ => {
                assert!(false, "incorrect error type");
            }
        }
    }

    /*
	#[test]
	fn test_pathology_flatlining() {

		// f(x)=1/e^(x^100) - 0.5, roots approx -0.996342 and 0.996342
		let f = |x: f64| x.powi(-100).exp() - 0.5;

		// d/dx = -100e^(-x^100) * x^99
		let df = |x: f64| -100.0 * (-x.powi(100)).exp() * x.powi(99);

        let root = newton_raphson(&f, &df, 0.99999, 100).expect("found root");
        assert!((root - 0.996342).abs() < 1e-6);
	}
*/
}
