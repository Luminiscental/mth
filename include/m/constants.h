#ifndef __m_constants_h__
#define __m_constants_h__

#include <cmath>
#include <limits>
#include <complex>

/* NOTE: Optional preprocessor flags:
 *
 *      M_ROW_MAJOR - Matrix values are stored row-major rather than column-major.
 *      M_ELIMINATION - Matrix inverses are calculated by Gaussian row elimination rather than by calculating the adjoint matrix.
 *      M_PRECISION - Value passed to std::setprecision (2 if not set).
 */

namespace m {

    template <typename T, size_t N>
    struct tvec; 

    constexpr std::complex<double> i = std::complex<double>(0, 1);

    template <typename T>
    constexpr T EPSILON = std::numeric_limits<T>::epsilon();

    template <typename T>
    constexpr T PI = static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640629);

    template <typename T>
    constexpr T TAU = static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257);

    template <typename T>
    constexpr tvec<T, 3> X_AXIS = tvec<T, 3>(1, 0, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Y_AXIS = tvec<T, 3>(0, 1, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Z_AXIS = tvec<T, 3>(0, 0, 1);

    // TODO: Move this somewhere else; "constants" doesn't really accurately describe it

    namespace util {

        template <typename T>
        bool checkZero(T x) {

            auto magnitude = std::abs(x);
            return magnitude <= EPSILON<decltype(magnitude)>;
        }

        template <typename T>
        bool checkEqual(T a, T b) {

            return checkZero(a - b);
        }
    }
}

#endif
