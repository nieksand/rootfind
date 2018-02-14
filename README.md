[![Build Status](https://travis-ci.org/nieksand/rootfind.svg?branch=master)](https://travis-ci.org/nieksand/rootfind) 
[![crates.io](https://img.shields.io/crates/v/rootfind.svg)](https://crates.io/crates/rootfind)
[![Released API docs](https://docs.rs/rootfind/badge.svg)](http://docs.rs/rootfind)

# Root Finding
***Work in progress.  Not ready for production use!***

Root finding algorithms implemented in Rust.

This package aims to provide robust numerical methods suitable for production
use.  It includes extensive documentation and test coverage.

Currently features:

* Bracket generation
* Bisection
* False Position, Illinois method

Some additional methods are only available in their "naive" form at this time.
These are suitable for reproducing results from academic literature but not for
production use:

* Newton-Raphson
* Halley's Method

Work is in progress on production-suitable variants which hybridize these
higher order methods with bisection to ensure convergence.

Custom convergence criteria can be supplied by the IsConverged trait.  Some
reasonable canned implementations are provided.

As with most numerical methods, root finding algorithms require that you
understand what you're trying to achieve, the nature of the input function, the
properties of the algorithm being used, and more.

Feedback is greatly appreciated.

# Usage
See the rustdocs for detailed documentation.

This quick example is an excerpt from tests/integration.rs.

    extern crate rootfind;
    
    use rootfind::bracket::{Bounds, BracketGenerator};
    use rootfind::solver::bisection;
    use rootfind::wrap::RealFn;
    
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

# Remaining Work
In terms of algorithms, there are three major things missing:

1. "Safe" variants of Newton-Raphson and Halley's Method which hybridize with a
   bracketing method to ensure global convergence.
2. Brent-Decker implementation for when no analytic derivatives are available.
3. Specialized routines for solving roots of Polynomials.

In terms of design, remaining work includes:

1. Allowing visibility into the solver state as it runs.
2. Potentially allowing optimized Newton-Raphson where the fraction f(x)/f'(x)
   is supplied directly rather than being computed at runtime.  Cancellation of
   terms provides an opportunity for performance optimization.

There are also two projects I want to cross-validate both implementations and
overall design against.  Specifically the C++ Boost Root Finders and the ones
supplied in the Gnu Scientific Library.  The Numerical Recipes book also has
solid implementations, but I want to avoid copyright issues so I'm mostly
staying away from it.

As expected, this project uses semantic versioning (major.minor.patch).  The
remaining work mostly falls under 'minor' increments.  When that's all done, I
would like some external review or feedback before cutting the official 1.0.0
release.

# References
The Numerical Recipes book covers both implementation and methodology for
root-finding in depth:

William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P.
Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific Computing
(3 ed.). Cambridge University Press, New York, NY, USA.

This is a top resource for practioners.  However, the code examples are
encumbered by copyright so the Rust rootfind library steers clear of NR's
implementations.

Another reasonable introduction to root finding can be found in:

Recktenwald, G. W. (2000). Numerical methods with MATLAB: implementations and
applications. Upper Saddle River, NJ: Prentice Hall.

Wikipedia's "Root-finding algorithm" page provides a high-level overview of
root-finding techniques, but it lacks the guidance and detail for practioners.
The algorithm specific pages are worth looking at.

I have also found the Boost, SciPy, and Gnu Scientific Library root-finding
implementations and documentation to be helpful.

# Author
This was written by Niek Sanders (niek.sanders@gmail.com).

# Unlicense
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
