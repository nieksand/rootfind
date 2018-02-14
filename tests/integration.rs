extern crate rootfind;

use rootfind::bracket::{Bounds, BracketGenerator};
use rootfind::solver::bisection;
use rootfind::wrap::RealFn;

#[test]
fn test_end_to_end() {
    // roots at 0, pi, 2pi, ...
    let f_inner = |x: f64| x.sin();

    // rootfind determines via traits what is f(x), df(x), d2f(x), etc.
    // the RealFn wrapper annotates our closure accordingly.
    let f = RealFn::new(&f_inner);

    // search for root-holding brackets
    let window_size = 0.1;
    let bounds = Bounds::new(-0.1, 6.3);

    for (i, b) in BracketGenerator::new(&f, bounds, window_size)
        .into_iter()
        .enumerate()
    {
        // find root using bisection method
        let max_iterations = 100;
        let computed_root = bisection(&f, &b, max_iterations).expect("found root");

        // demonstrate that we found root
        let pi = std::f64::consts::PI;
        let expected_root = (i as f64) * pi;

        assert!(
            (computed_root - expected_root).abs() < 1e-9,
            format!("got={}, wanted={}", computed_root, expected_root)
        );
    }
}
