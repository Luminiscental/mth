#ifndef __m_mat_h__
#define __m_mat_h__

#include <iostream>
#include <array>
#include <type_traits>

#include <m/vec.h>
#include <m/quat.h>
#include <m/comp.h>

// NOTE: Safeguard so mat_content.h doesn't get included anywhere else
#define __m_mat_content_toggle__ 

namespace m {

    template <typename T, size_t N, size_t M>
    class tmat;

    template <typename T, size_t N, size_t M>
    tmat<T, N, M> operator+(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    tmat<T, N, M> operator-(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    tmat<T, N, M> operator*(const T &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    tmat<T, N, M> operator*(const tmat<T, N, M> &lhs, const T &rhs);

    template <typename T, size_t N, size_t M>
    tmat<T, N, M> operator/(const tmat<T, N, M> &lhs, const T &rhs);

    template <typename T, size_t N, size_t M, size_t O>
    tmat<T, O, M> operator*(const tmat<T, N, M> &lhs, const tmat<T, O, N> &rhs);

    template <typename T, size_t N>
    tmat<T, N, N> operator/(const tmat<T, N, N> &lhs, const tmat<T, N, N> &rhs);

    template <typename T, size_t N, size_t M>
    tvec<T, M> operator*(const tmat<T, N, M> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N, size_t M>
    bool operator==(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    bool operator!=(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    // TODO: Either make this less crappy or remove it
    template <typename T, size_t N, size_t M>
    std::ostream &operator<<(std::ostream &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, typename A>
    class tmat_aug;

    template <typename T, size_t N, typename A>
    std::ostream &operator<<(std::ostream &lhs, const tmat_aug<T, N, A> &rhs);

    template <typename T, size_t N, typename A>
    class tmat_aug {

    private:

        std::array<A, N> aux;

        // TODO: Augment non-square matrices?
        tmat<T, N, N> matrix;

        void eliminateFromBelow(size_t x, size_t y);
        void eliminateFromRight(size_t x, size_t y);

    public:

        tmat_aug() noexcept;
        tmat_aug(const tmat<T, N, N> &matrix, const std::array<A, N> &aux) noexcept;

        tmat<T, N, N> coefficients() const; 
        std::array<A, N> auxilary() const;

        // TODO: degenerate solution cases
        std::array<A, N> solve() const;

        // NOTE: returns N on empty rows
        size_t leadingIndex(size_t row) const;
        // NOTE: returns zero on empty rows
        T leadingValue(size_t row) const;

        bool columnIsZero(size_t x) const;
        bool rowIsZero(size_t y) const;
        bool hasZeroRow() const;
        bool singular() const;

        void swapRows(size_t a, size_t b); 
        void scaleRow(size_t index, T scalar);
        void addRow(size_t targetRow, size_t sourceRow, T scalar = 1); 

        void setRow(size_t index, const tvec<T, N> &val, A auxVal); 

        tmat_aug<T, N, A> ordered() const;
        tmat_aug<T, N, A> rowEchelon() const; 
        tmat_aug<T, N, A> reducedRowEchelon() const; 

        // TODO: Arithmetic on augmented matrices?
    };   

#define M 1
#define N 1

    // NOTE: 1,1
    template <typename T>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N
#undef M

#define N 1

    // NOTE: 1,M
    template <typename T, size_t M>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N

#define M 1

    // NOTE: N,1
    template <typename T, size_t N>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef M

    // NOTE: N,M
    template <typename T, size_t N, size_t M>
    class tmat {

#include <m/mat_content.h>

    };

    // NOTE: Specialized static functions 
    namespace mat {

        template <typename T>
        tmat<T, 4, 4> scale(const tvec<T, 3> &factors);

        template <typename T>
        tmat<T, 4, 4> scale(T factor);

        template <typename T>
        tmat<T, 4, 4> translation(const tvec<T, 3> &offset);

        template <typename T>
        tmat<T, 4, 4> rotation(const tquat<T> &rep);

        template <typename T>
        tmat<T, 4, 4> rotation(T angle, const tvec<T, 3> &axis);

        // TODO: Ortho and perspective projections
    }

#define CREATE_SQR_ALIASES(n)       using imat ## n = tmat<int, n, n>; \
                                    using lmat ## n = tmat<long, n, n>; \
                                    using  mat ## n = tmat<float, n, n>; \
                                    using dmat ## n = tmat<double, n, n>; \
                                    using cmat ## n = tmat<comp, n, n>;

#define CREATE_RECT_ALIASES(n, m)   using imat ## n ## x ## m = tmat<int, n, m>; \
                                    using lmat ## n ## x ## m = tmat<long, n, m>; \
                                    using  mat ## n ## x ## m = tmat<float, n, m>; \
                                    using dmat ## n ## x ## m = tmat<double, n, m>; \
                                    using cmat ## n ## x ## m = tmat<comp, n, m>;

#define CREATE_ALIASES(n)   CREATE_SQR_ALIASES(n) \
                            CREATE_RECT_ALIASES(n, 1) \
                            CREATE_RECT_ALIASES(n, 2) \
                            CREATE_RECT_ALIASES(n, 3) \
                            CREATE_RECT_ALIASES(n, 4) \
                            CREATE_RECT_ALIASES(n, 5) \
                            CREATE_RECT_ALIASES(n, 6) \
                            CREATE_RECT_ALIASES(n, 7) \
                            CREATE_RECT_ALIASES(n, 8) \
                            CREATE_RECT_ALIASES(n, 9)

    CREATE_ALIASES(1)
    CREATE_ALIASES(2)
    CREATE_ALIASES(3)
    CREATE_ALIASES(4)
    CREATE_ALIASES(5)
    CREATE_ALIASES(6)
    CREATE_ALIASES(7)
    CREATE_ALIASES(8)
    CREATE_ALIASES(9)

#undef CREATE_ALIASES
#undef CREATE_SQR_ALIASES

}

#define N 1
#define M 1

// NOTE: 1,1
#include <m/mat_impl.h>

#undef M
#undef N

#define N 1

// NOTE: 1,M
#include <m/mat_impl.h>

#undef N

#define M 1

// NOTE: N,1
#include <m/mat_impl.h>

#undef M

// NOTE: N,M
#include <m/mat_impl.h>

#undef __m_mat_content_toggle__

#endif
