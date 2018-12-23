#ifndef __Maths_mat_h__
#define __Maths_mat_h__

#include <string>
#include <complex>
#include <memory>
#include <iostream>
#include <iomanip>

#ifndef __Maths_vec_h__
#include <Maths/vec.h>
#endif

#ifndef __Maths_constants_h__
#include <Maths/constants.h>
#endif

#ifndef __Maths_quat_h__
#include <Maths/quat.h>
#endif

#define __Maths_mat_content_toggle__ // Safeguard so mat_content.h doesn't get included anywhere else

namespace m {

    // TODO: Non-square matrices
    
    template <typename T, size_t N>
    struct tmat;

    // tmat<T, 2>

#define __Maths_mat_basecaseimpl__

#define N 2

    template <typename T>
    struct tmat<T, 2> {

#include <Maths/mat_content.h>

    };

#undef N

#undef __Maths_mat_basecaseimpl__

    // tmat<T, N> for N > 2

    template <typename T, size_t N>
    struct tmat {

#include <Maths/mat_content.h>

    };

    // transformations
    
    namespace mat {

        template <typename T>
        tmat<T, 4> scale(const tvec<T, 3> &factors) {

            return tmat<T, 4>{factors.get(0), 0,              0,              0,
                              0,              factors.get(1), 0,              0,
                              0,              0,              factors.get(2), 0,
                              0,              0,              0,              1};
        }

        template <typename T>
        tmat<T, 4> scale(T factor) {

            return scale(tvec<T, 3>{factor, factor, factor});
        }

        template <typename T>
        tmat<T, 4> translate(const tvec<T, 3> &offset) {

            return tmat<T, 4>{1, 0, 0, 0, // NOTE: Layout transposed because column-major
                              0, 1, 0, 0,
                              0, 0, 1, 0,
                              offset.x(), offset.y(), offset.z(), 1};
        }

        template <typename T>
        tmat<T, 4> rotate(T angle, const tvec<T, 3> &axis) {

            return rotate(tquat<T>::rotation(angle, axis));
        }

        template <typename T>
        tmat<T, 4> rotate(const tquat<T> &rep) {
            
            tvec<T, 3> rotatedX = rep.rotate(m::x_axis<T>);
            tvec<T, 3> rotatedY = rep.rotate(m::y_axis<T>);
            tvec<T, 3> rotatedZ = rep.rotate(m::z_axis<T>);

            return tmat<T, 4>{rotatedX.x(), rotatedX.y(), rotatedX.z(), 0, // NOTE: Layout transposed because column-major
                              rotatedY.x(), rotatedY.y(), rotatedY.z(), 0,
                              rotatedZ.x(), rotatedZ.y(), rotatedZ.z(), 0,
                              0,            0,            0,            1};
        }

        // TODO: Ortho and perspective projections
    }

#define TYPEDEF_MAT(n) typedef tmat<int, n> imat ## n; \
                       typedef tmat<long, n> lmat ## n; \
                       typedef tmat<float, n> mat ## n; \
                       typedef tmat<double, n> dmat ## n; \
                       typedef tmat<std::complex<double>, n> cmat ## n;

    TYPEDEF_MAT(2)
    TYPEDEF_MAT(3)
    TYPEDEF_MAT(4)

#undef TYPEDEF_MAT

}

#undef __Maths_mat_content_toggle__

#endif
