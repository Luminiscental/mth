
# m - A Maths library

No one wanted it; no one needed it, but it's here now.

m is a maths library built around template classes to work with objects such as vectors, matrices and quaternions, with helper classes to deal with solving equations. No guarantees of stability or quality atm so be warned.

For example code see src/main.cpp

### Optional preprocessor flags:

* `m_ROW_MAJOR` - When set the values of matrices are stored row major rather than column-major. This doesn't affect behaviour when using relevant functions / classes, only the implementation.
* `m_ELIMINATION` - When set the `inverse` method of `m::tmat` uses Gaussian elimination instead of recursing on minor matrices to find the adjoint.
* `m_PRECISION` - If set this is passed to `std::setprecision` in pretty-print functions; it's the number of decimal places displayed.

### Features:

* Relevant constants defined: e.g. pi and tau, coordinate bases.
* N-dimensional vector of any given scalar type (that has most arithmetic operators / functions defined on it) for N at least 1.
* N,M-dimensional matrix of any given scalar type (same as above) for N and M at least 2.
* Augmentations of square matrices with any type that supports addition and scalar multiplication (e.g. `float`, `m::comp`, `m::tvec`), with row operations and functions to convert to echelon / reduced echelon form.
* Quaternions of any given scalar type (you get the idea) with arithmetic operators and functions for representing rotations in 3-space (and converting to the equivalent matrices).
* Complex number representations with arithmetic and overloads for `std::exp`, `std::cos`, `std::sin` and `std::abs`.
* Polynomials with arithmetic, a `solve()` method for finding roots and a `value()` method for evaluating.
* Differentiation and integration of polynomials.
* Somewhat pretty printing.

### Planned features:

* Numeric methods for finding roots of arbitrary functions
* Numeric differentiation and integration of arbitrary complex functions.
* Numerical series expansions
