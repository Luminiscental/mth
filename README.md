
# m - A Maths library

No one wanted it; no one needed it, but it's here now.

m is a maths library built around template classes to work with objects such as vectors, matrices and quaternions, with helper classes to deal with solving equations. Currently a header library because of templates but that may not necessarily be the case in the future.

### Optional preprocessor flags:

* `M_ROW_MAJOR` - When set the values of matrices are stored row major rather than column-major. This isn't relevant for specifying values in constructors as they are always interpreted as row-major but it means that a list of vectors is interpreted as rows rather than columns.
* `M_ELIMINATION` - When set the `inverse` method of `tmat` uses Gaussian elimination instead of recursing on minor matrices to find the adjoint.
* `M_PRECISION` - If set this is passed to `std::setprecision` in pretty-print functions; it's the number of decimal places displayed.
* More - soon tm

### Features:

* Relevant constants defined: e.g. pi and tau, coordinate bases.
* N-dimensional vector of any given scalar type (that has most arithmetic operators / functions defined on it) for N = 1, 2, 3, ...
* N-dimensional square matrix of any given scalar type (same as above) for N = 2, 3, 4, ...
* Augmentations of square matrices with any type that supports addition and scalar multiplication (e.g. `float`, `std::complex`, `tvec`), with row operations and functions to convert to echelon / reduced echelon form.
* Quaternions of any given scalar type (you get the idea) with arithmetic operators and functions for representing rotations in 3-space (and converting to the equivalent matrices).
* Classes for solving certain types of equations, e.g. linear systems, polynomials for given scalar types (currently incomplete / WIP).
* Somewhat pretty printing.

### Planned features:

* Complex numbers because redundancy is my middle name
* Numeric methods for finding roots of arbitrary functions
* Analytic calculus on polynomials and numeric on arbitrary functions
* Numerical series expansions

### Example code:

```c

#define M_PRECISION 3

#include <iostream>

#include <m/constants.h>
#include <m/vec.h>
#include <m/mat.h>
#include <m/quat.h>

int main() {

    m::mat4 transformation = m::mat::rotate(m::PI<float> / 3, m::X_AXIS<float>) * m::mat::translate(m::vec3(1.0f, 13.0f, -2.0f));

    m::vec3 initialPosition(1, 2, 1);
    m::vec4 positionHandle(initialPosition.x(), initialPosition.y(), initialPosition.z(), 1);
    
    m::vec3 transformedPosition = (transformation * positionHandle).xyz();

    m::quat additionalRotation = m::quat::rotation(m::TAU<float> / 12, m::vec3(0.5f, 0.5f, -0.1f));

    additionalRotation *= m::quat(2, 3, 1, 5);
    additionalRotation.real() -= 4;

    transformedPosition = additionalRotation.rotate(transformedPosition);

    std::cout << std::endl << "Transformation matrix:" << std::endl;
    std::cout << std::endl << m::mat::rotate(additionalRotation) * transformation << std::endl;

    std::cout << std::endl << "Went from " << initialPosition << " to " << transformedPosition << std::endl;

    return 0;
}

/* Output:

    Transformation matrix:

    |   -0.784  -0.557  0.190   0.000   |
    |   0.027   -0.207  0.416   0.000   |
    |   -0.442  0.255   -0.702  0.000   |
    |   1.000   13.000  -2.000  1.000   |

    Went from (1.000, 2.000, 1.000) to (-0.657, -0.073, -0.634)

 */
```
