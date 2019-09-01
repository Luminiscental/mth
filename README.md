
# mth - A Maths library

mth is a maths library built around template classes to work with objects such as vectors, matrices
and quaternions, with helper classes to deal with solving equations. No guarantees of stability or
quality atm so be warned.

### Building:

To build from source, you will need cmake and a c++ compiler supporting c++17.

First, clone the repo and run cmake:
```
$ git clone https://github.com/Luminiscental/mth
$ cd mth && mkdir build && cd build
$ cmake .. # Optionally pass -DTESTS=ON to run tests when building
```

Then run the generated build files, on mac/linux there should be a Makefile in the build directory
that you can run with `make`. On windows there should be a Visual Studio solution file `mth.sln`,
open this and build the solution.

After building there should be a `libmth.a` or `mth.lib` file for linking into programs that use
mth.

### Features:

* Relevant constants defined: e.g. pi and tau, coordinate bases.
* N-dimensional vector of any given scalar type (that has most arithmetic operators / functions
  defined on it) for N at least 1.
* N,M-dimensional matrix of any given scalar type (same as above) for N and M at least 2.
* Augmentations of square matrices with any type that supports addition and scalar multiplication
  (e.g. `float`, `mth::comp`, `mth::tvec`), with row operations and functions to convert to echelon
  / reduced echelon form.
* Quaternions of any given scalar type (you get the idea) with arithmetic operators and functions
  for representing rotations in 3-space (and converting to the equivalent matrices).
* Complex number representations with arithmetic and overloads for `std::exp`, `std::cos`,
  `std::sin` and `std::abs`.
* Polynomials with arithmetic, a `solve()` method for finding roots and a `value()` method for
  evaluating.
* Numerical calculation of the limits of sequences or of complex functions at a point.
* Power series with complex coefficients allowing evaluation at points and
  differentiation/integration.
* Differentiation and integration of polynomials.
* Somewhat pretty printing.

### Planned features:

* Numeric methods for finding roots of arbitrary functions
* Numerical series expansions
