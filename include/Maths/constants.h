#ifndef __Maths_constants_h__
#define __Maths_constants_h__

#include <cmath>
#include <limits>

#ifndef __Maths_vec_h__
#include <Maths/vec.h>
#endif

#define Maths_EPSILON(T) std::numeric_limits<T>::epsilon()

#define Maths_PI(T)  static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640629)
#define Maths_TAU(T) static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257)

namespace m {

    template <typename T, size_t N>
    struct tvec; // idk why this is needed

    std::complex<double> i(0, 1);

    template <typename T>
    tvec<T, 3> x_axis{1, 0, 0};
                                
    template <typename T>
    tvec<T, 3> y_axis{0, 1, 0};
                                
    template <typename T>
    tvec<T, 3> z_axis{0, 0, 1};

    namespace util {

        template <typename T>
        bool checkZero(T x) {

            return std::abs(x) < Maths_EPSILON(T);
        }
    }
}

#endif
