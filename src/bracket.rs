use std::f64;

/// Bounds represents the closed interval [a,b].
#[derive(Clone, Debug, PartialEq)]
pub struct Bounds {
    pub a: f64,
    pub b: f64,
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
pub fn is_sign_change(lhs: f64, rhs: f64) -> bool {
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
        Bounds::new(f64::NAN, -2.0);
    }

    #[test]
    #[should_panic]
    fn test_bounds_new_infinite() {
        Bounds::new(f64::NEG_INFINITY, f64::INFINITY);
    }

    #[test]
    fn test_bracket_generator_hits() {
        let f = |x: f64| x.sin();
        let pi = f64::consts::PI;
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

    #[test]
    fn test_first_bracket_even_degree() {
        // root at zero is of even degree (touches but does not cross x-axis).
        // bracketing won't find that.
        let f = |x| x * x;
        let win = first_bracket(&f, &Bounds::new(-4.5, 4.5), 1.0);
        assert!(win.is_none());
    }
}
