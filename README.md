
# m - A Maths library

No one wanted it; no one needed it, but it's here now.

m is a maths library built around template classes to work with objects such as vectors, matrices and quaternions, with helper classes to deal with solving equations. Currently a header library because of templates but that may not necessarily be the case in the future.

### Optional preprocessor flags:

* `m_ROW_MAJOR` - When set the values of matrices are stored row major rather than column-major. This isn't relevant for specifying values in constructors as they are always interpreted as row-major but it means that a list of vectors is interpreted as rows rather than columns.
* `m_ELIMINATION` - When set the `inverse` method of `tmat` uses Gaussian elimination instead of recursing on minor matrices to find the adjoint.
* `m_PRECISION` - If set this is passed to `std::setprecision` in pretty-print functions; it's the number of decimal places displayed.
* More - soon tm

### Features:

* Relevant constants defined: e.g. pi and tau, coordinate bases.
* N-dimensional vector of any given scalar type (that has most arithmetic operators / functions defined on it) for N at least 1.
* N,M-dimensional matrix of any given scalar type (same as above) for N and M at least 2.
* Augmentations of square matrices with any type that supports addition and scalar multiplication (e.g. `float`, `std::complex`, `tvec`), with row operations and functions to convert to echelon / reduced echelon form.
* Quaternions of any given scalar type (you get the idea) with arithmetic operators and functions for representing rotations in 3-space (and converting to the equivalent matrices).
* Polynomials with a `solve()` method for finding roots and a `value()` method for evaluating.
* Somewhat pretty printing.

### Planned features:

* Complex numbers because redundancy is my middle name
* Numeric methods for finding roots of arbitrary functions
* Analytic calculus on polynomials and numeric on arbitrary functions
* Numerical series expansions

### Example code:

main.cpp:

```c

#include <iostream>

#include <m/m.h>

#include <m/vec.h>
#include <m/mat.h>
#include <m/quat.h>

int main() {

    m::mat4 transformation = m::mat::rotation(m::PI<float> / 3, m::X_AXIS<float>)
                           * m::mat::translate(m::vec3(1.0f, 13.0f, -2.0f));

    m::vec3 initialPosition(1, 2, 1);
    m::vec4 positionHandle(initialPosition.x(), initialPosition.y(), initialPosition.z(), 1);
    
    m::vec3 transformedPosition = (transformation * positionHandle).xyz();

    m::quat additionalRotation = m::quat::rotation(m::TAU<float> / 12, m::vec3(0.5f, 0.5f, -0.1f));

    additionalRotation *= m::quat(2, 3, 1, 5);
    additionalRotation.real() -= 4;

    transformedPosition = additionalRotation.rotate(transformedPosition);

    std::cout << std::endl << "Transformation matrix:" << std::endl;
    std::cout << std::endl << m::mat::rotation(additionalRotation) * transformation << std::endl;

    std::cout << std::endl << "Went from " << initialPosition << " to " << transformedPosition << std::endl;

    m::tmat<double, 8, 8> bigMatrix(0.0, 3.0, -1.0, 5.0, 7.0, 6.0, 2.0, -3.0,
                                    1.0, 5.0, 2.0, -3.0, 2.0, -2.0, -2.0, 5.0,
                                    7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,
                                    3.0, 2.0, -1.0, 3.0, -5.0, -6.0, -1.0, 1.0,
                                    5.0, 3.0, 4.0, 3.0, 5.0, -2.0, -3.0, 1.0,
                                    6.0, 4.0, 3.0, 2.0, 7.0, -5.0, -7.0, 11.0,
                                    0.0, 2.0, 0.0, 1.5, 2.5, -3.0, 5.0, 4.0,
                                    2.0, 4.0, 8.0, 12.0, -11.0, 2.0, 5.0, -6.0);

    std::cout << std::endl << "A : " << std::endl;
    std::cout << std::endl << bigMatrix << std::endl;

    std::cout << std::endl << "A^{-1} : " << std::endl;
    std::cout << std::endl << bigMatrix.inverse() << std::endl;

    std::cout << std::endl;

    return 0;
}

```

Output:

```

Transformation matrix:

|	-0.78	-0.21	0.52	-4.53	|
|	-0.32	-0.21	-0.16	-2.69	|
|	-0.27	-0.51	-0.70	-5.50	|
|	0.00	0.00	0.00	1.00	|

Went from (1.00, 2.00, 1.00) to (-5.21, -3.58, -7.49)

A : 

|	0.00	3.00	-1.00	5.00	7.00	6.00	2.00	-3.00	|
|	1.00	5.00	2.00	-3.00	2.00	-2.00	-2.00	5.00	|
|	7.00	6.00	5.00	4.00	3.00	2.00	1.00	0.00	|
|	3.00	2.00	-1.00	3.00	-5.00	-6.00	-1.00	1.00	|
|	5.00	3.00	4.00	3.00	5.00	-2.00	-3.00	1.00	|
|	6.00	4.00	3.00	2.00	7.00	-5.00	-7.00	11.00	|
|	0.00	2.00	0.00	1.50	2.50	-3.00	5.00	4.00	|
|	2.00	4.00	8.00	12.00	-11.00	2.00	5.00	-6.00	|

A^{-1} : 

|	-0.07	-0.13	0.20	0.05	-0.08	0.02	-0.02	-0.06	|
|	0.09	0.18	-0.02	0.08	-0.00	-0.05	-0.03	0.01	|
|	-0.08	0.03	-0.04	-0.13	0.13	-0.02	0.03	0.06	|
|	0.07	-0.06	-0.07	0.02	-0.00	0.07	0.01	0.05	|
|	0.03	-0.00	-0.04	-0.03	0.12	-0.03	0.05	-0.02	|
|	0.02	-0.03	0.10	-0.06	-0.19	0.07	-0.07	0.01	|
|	-0.04	-0.04	0.08	-0.02	-0.03	-0.04	0.12	-0.01	|
|	-0.02	-0.03	0.04	-0.05	-0.17	0.12	0.03	0.02	|

```
