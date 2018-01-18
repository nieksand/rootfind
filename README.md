[![Build Status](https://travis-ci.org/nieksand/rootfind.svg?branch=master)](https://travis-ci.org/nieksand/rootfind)

# Root Finding

***Work in progress.  Not ready for use!***

Root finding functions implemented in Rust.

Currently features:

* Newton-Raphson
* Halley's Method
* Bisection
* Bracket generation

# General Approach

1. Know what you are trying to achieve.
   1. How much accuracy do you need in your final answer?
   2. Do you need to trade off accuracy for speed?
   3. Are there any particular roots your interested in?
   4. Do you care about complex-valued roots?
   5. Is there an existing, analytic solution you can use over a numerical one?

2. Understand the function or family of functions you are trying to solve.
   1. Plot your function or some samples from the family!
   1. Are they continuous in the region of interest?
   2. Are there any singularities in the region of interest?
   3. Is there a possibility of distinct yet closely-spaced roots?  
      Do you care if they're treated as just one?

...

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
