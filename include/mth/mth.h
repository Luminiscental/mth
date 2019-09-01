#ifndef mth_mth_h__
#define mth_mth_h__

/* <mth/mth.h> - main header
 *      This includes general constant definitions and utility functions.
 */

// Library info

/* File naming schemes:
 *      name.h - standard header file with class / function declarations
 *          (manually include these)
 *      name_impl.h - header file implementing template functions / classes
 *          (included already in name.h)
 *      name_content.h - header file containing the content of a class to be
 *          included multiple times by header.h for template specialization
 *          (included already in name.h)
 *      name.cpp - source file implementing non-template functions / classes
 *          (compiled into libM)
 */

/* Class naming schemes:
 *      tname - Basic data / value storage class templated for a scalar type T
 *      Name - More advanced structure / object; non-template class
 */

/* Optional preprocessor flags:
 *
 *      mth_ROW_MAJOR - Matrix values are stored row-major rather than
 *          column-major.
 *      mth_ELIMINATION - Matrix inverses are calculated by Gaussian row
 *          elimination using tmat_aug rather than by calculating the
 *          adjoint matrix.
 */

#include <cmath>
#include <limits>

// TODO: Sequence wrapper class (e.g. for recursive sequences)
// TODO: Rational functions
// TODO: Laurent series
// TODO: Better exception / debug info maybe

namespace mth {

    // Smallest value above zero; used for equality conditions on non-exact types
    template <typename T>
    constexpr T epsilon = std::numeric_limits<T>::epsilon();

    // Mathematical scalar constants

    template <typename T>
    constexpr T pi =  static_cast<T>(3.141592653589793238462643383279502884197169399375100582097494459230781640629);
    
    template <typename T>
    constexpr T tau = static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257);

    template <typename T>
    constexpr T e =   static_cast<T>(2.71828182845904523536028747135266249775724709369995957496696762772407663035);

    double factorial(size_t n);

    // Global utility functions

    namespace util {

        // Check whether a value is equal to zero
        template <typename T>
        constexpr bool isZero(const T &x) noexcept {

            return x == static_cast<T>(0);
        }

        // For floating points check as an inequality with the smallest value
        
        template <>
        constexpr bool isZero(const float &x) noexcept {

            // for now cmath functions aren't constexpr
            /*
            using std::abs;
            auto magnitude = abs(x);
            */
            auto magnitude = x < 0 ? -x : x;
            return magnitude <= mth::epsilon<decltype(magnitude)>;
        }

        template <>
        constexpr bool isZero(const double &x) noexcept {

            // for now cmath functions aren't constexpr
            /*
            using std::abs;
            auto magnitude = abs(x);
            */
            auto magnitude = x < 0 ? -x : x;
            return magnitude <= mth::epsilon<decltype(magnitude)>;
        }

        // Check whether two values are equal
        template <typename T>
        constexpr bool isEqual(const T &a, const T &b) noexcept {

            return a == b;
        }

        // For floating points defer to isZero on a - b

        template <>
        constexpr bool isEqual(const double &a, const double &b) noexcept {

            return isZero(a - b);
        }
        
        template <>
        constexpr bool isEqual(const float &a, const float &b) noexcept {

            return isZero(a - b);
        }
    }
}

#endif
