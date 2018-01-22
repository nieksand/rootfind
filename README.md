[![Build Status](https://travis-ci.org/nieksand/rootfind.svg?branch=master)](https://travis-ci.org/nieksand/rootfind) 
[![crates.io](https://img.shields.io/crates/v/rootfind.svg)](https://crates.io/crates/rootfind)

# Root Finding
***Work in progress.  Not ready for use!***

Root finding algorithms implemented in Rust.

Currently features:

* Newton-Raphson
* Halley's Method
* Bisection
* Bracket generation

This package aims to provide robust algorithms for finding real roots of
univariate functions.

Custom convergence criteria can be supplied by the IsConverged trait.  Some
reasonable canned implementations are provided.

As with most numerical methods, root finding algorithms require that you
understand what you're trying to achieve, the nature of the input function, the
properties of the algorithm being used, and more.

The wikipedia page on "Root-finding algorithm" is a reasonable introduction.

Feedback is greatly appreciated.

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

Then in the 'misc' bucket of work:

1. Flesh out README document to include working examples.
2. Flesh out rustdocs and maintain a published version, linked via badge.
3. Add crates.io version badge.

As expected, this project uses semantic versioning (major.minor.patch).  The
remaining work mostly falls under 'minor' increments.  When that's all done, I
would like some external review or feedback before cutting the official 1.0.0
release.

# Usage
...

# References
A reasonable introduction to root finding can be found in:

Recktenwald, G. W. (2000). Numerical methods with MATLAB: implementations and
applications. Upper Saddle River, NJ: Prentice Hall.

The Wikipedia page on "Root-finding algorithm" is also a good overview.  The
quality of various algorithm-specific pages varies dramatically.

I have also found both the C++ Boost and Gnu Scientific Library root-finding
implementations and documentation to be quite helpful.

The Numerical Recipes book covers implementation in detail and is available to
read online, but the implementations are encumbered by copyright.

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
