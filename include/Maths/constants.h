#ifndef __Maths_constants_h__
#define __Maths_constants_h__

#include <limits>

#define Maths_EPSILON(T) std::numeric_limits<T>::epsilon()

#define Maths_PI(T)  static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640629)
#define Maths_TAU(T) static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257)

namespace m {

    namespace util {

        // requires std::abs(T), std::numeric_limits<T>::epsilon()
        template <typename T>
        bool checkZero(T x);
    }
}

// Template implementation

template <typename T>
bool m::util::checkZero(T x) {

    return std::abs(x) < Maths_EPSILON(T);
}

#endif
