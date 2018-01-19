[![Build Status](https://travis-ci.org/nieksand/rootfind.svg?branch=master)](https://travis-ci.org/nieksand/rootfind)

# Root Finding
***Work in progress.  Not ready for use!***

Root finding functions implemented in Rust.

Currently features:

* Newton-Raphson
* Halley's Method
* Bisection
* Bracket generation

This package aims to provide robust algorithms for finding real roots of
univariate functions.

As with most numerical methods, root finding algorithms require that you
understand your goals, the nature of the input function, the properties of the
algorithm being used, and more.

Feedback is greatly appreciated.

# Usage
...

# References
A good introduction for Newton-Raphson is in:

    Richard W Hamming. 1986. Numerical Methods for Scientists and Engineers 
    (2nd Ed.). Dover Publications, Inc., New York, NY, USA.

The "Root-finding algorithm" wikipedia page is good topic intro.

The "Newton's method" wikipedia page is also quite comprehensive.

The Numerical Recipes book also has a strong treatment, but it is encumbered by
copyright when it comes to implementation.

The boost algorithms page has a decent collection and overview:

    http://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots.html

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
