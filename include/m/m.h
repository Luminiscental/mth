#ifndef __m_m_h__
#define __m_m_h__

#define __m_impl__

/* NOTE: Optional preprocessor flags:
 *
 *      m_ROW_MAJOR - Matrix values are stored row-major rather than column-major.
 *      m_ELIMINATION - Matrix inverses are calculated by Gaussian row elimination rather than by calculating the adjoint matrix.
 *      m_PRECISION - Value passed to std::setprecision (2 if not set).
 */

#ifndef m_PRECISION
#define m_PRECISION 2
#endif

namespace m {

    template <typename T>
    constexpr T EPSILON;

    template <typename T>
    constexpr T PI;

    template <typename T>
    constexpr T TAU;

    namespace util {

        template <typename T>
        auto checkZero(const T &x);

        template <typename T>
        auto checkEqual(const T &a, const T &b);
    }
}

#include <m/m_impl.h>

#endif
