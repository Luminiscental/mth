#ifndef __Maths_constants_h__
#define __Maths_constants_h__

#include <limits>

#define Maths_EPSILON(T) std::numeric_limits<T>::epsilon()

namespace m {

    namespace util {

        template <typename T> // requires std::abs(T)
        bool checkZero(T x);
    }
}

// Template implementation

template <typename T>
bool m::util::checkZero(T x) {

    return std::abs(x) < Maths_EPSILON(T);
}

#endif
