#ifndef __Maths_mat_h__
#define __Maths_mat_h__

#include <string>
#include <complex>
#include <memory>

#ifndef __Maths_vec_h__
#include <Maths/vec.h>
#endif

#ifndef __Maths_constants_h__
#include <Maths/constants.h>
#endif

#define __Maths_mat_content_toggle__ // Safeguard so mat_content.h doesn't get included anywhere else

namespace m {

    template <typename T, size_t N>
    struct tmat;

#define __Maths_mat_basecaseimpl__

#define N 2

    template <typename T>
    struct tmat<T, 2> {

#include <Maths/mat_content.h>

    };

#undef N

#undef __Maths_mat_basecaseimpl__

    template <typename T, size_t N>
    struct tmat {

#include <Maths/mat_content.h>

    };

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
