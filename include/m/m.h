#ifndef __m_m_h__
#define __m_m_h__

/* NOTE: Optional preprocessor flags:
 *
 *      m_ROW_MAJOR - Matrix values are stored row-major rather than column-major.
 *      m_ELIMINATION - Matrix inverses are calculated by Gaussian row elimination rather than by calculating the adjoint matrix.
 *      m_PRECISION - Value passed to std::setprecision (2 if not set).
 */

#ifndef m_PRECISION
#define m_PRECISION 2
#endif

#include <limits>

namespace m {

    template <typename T>
    constexpr T EPSILON = std::numeric_limits<T>::epsilon();

    template <typename T>
    constexpr T PI = static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640629);
    
    template <typename T>
    constexpr T TAU = static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257);

    namespace util {

        template <typename T>
        bool checkZero(const T &x);

        template <typename T>
        bool checkEqual(const T &a, const T &b);
    }
}

#include <m/m_impl.h>

#endif
