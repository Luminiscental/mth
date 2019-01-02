#ifndef __m_mat_h__
#define __m_mat_h__

#include <iostream>
#include <array>
#include <type_traits>

#include <m/vec.h>
#include <m/quat.h>
#include <m/comp.h>

#define __m_mat_content_toggle__ // NOTE: Safeguard so mat_content.h doesn't get included anywhere else

namespace m {

    template <typename T, size_t N, size_t M>
    class tmat;

    template <typename T, size_t N, size_t M>
    auto operator+(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    auto operator-(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    auto operator*(const T &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, size_t M>
    auto operator*(const tmat<T, N, M> &lhs, const T &rhs);

    template <typename T, size_t N, size_t M>
    auto operator/(const tmat<T, N, M> &lhs, const T &rhs);

    template <typename T, size_t N, size_t M, size_t O>
    auto operator*(const tmat<T, N, M> &lhs, const tmat<T, O, N> &rhs);

    template <typename T, size_t N>
    auto operator/(const tmat<T, N, N> &lhs, const tmat<T, N, N> &rhs);

    template <typename T, size_t N, size_t M>
    auto operator*(const tmat<T, N, M> &lhs, const tvec<T, N> &rhs);

    // TODO: Either make this less crappy or remove it
    template <typename T, size_t N, size_t M>
    auto &operator<<(std::ostream &lhs, const tmat<T, N, M> &rhs);

    template <typename T, size_t N, typename A>
    class tmat_aug;

    template <typename T, size_t N, typename A>
    auto &operator<<(std::ostream &lhs, const tmat_aug<T, N, A> &rhs);

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

        auto coefficients() const; 
        auto auxilary() const;

        // TODO: degenerate solution cases
        auto solve() const;

        // NOTE: returns N on empty rows
        auto leadingIndex(size_t row) const;
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

        auto ordered() const;
        auto rowEchelon() const; 
        auto reducedRowEchelon() const; 
    };   

    // NOTE: 2,2

#define M 2
#define N 2

    template <typename T>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N
#undef M

    // NOTE: 2,M

#define N 2

    template <typename T, size_t M>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N

    // NOTE: N,2

#define M 2

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
        auto scale(const tvec<T, 3> &factors);

        template <typename T>
        auto scale(T factor);

        template <typename T>
        auto translate(const tvec<T, 3> &offset);

        template <typename T>
        auto rotation(const tquat<T> &rep);

        template <typename T>
        auto rotation(T angle, const tvec<T, 3> &axis);

        // TODO: Ortho and perspective projections
    }

#define CREATE_ALIASES(n) using imat ## n = tmat<int, n, n>; \
                          using lmat ## n = tmat<long, n, n>; \
                          using  mat ## n = tmat<float, n, n>; \
                          using dmat ## n = tmat<double, n, n>; \
                          using cmat ## n = tmat<m::comp, n, n>;

    CREATE_ALIASES(2)
    CREATE_ALIASES(3)
    CREATE_ALIASES(4)
    CREATE_ALIASES(5)
    CREATE_ALIASES(6)
    CREATE_ALIASES(7)
    CREATE_ALIASES(8)
    CREATE_ALIASES(9)

#undef CREATE_ALIASES

}

#define N 2
#define M 2

// NOTE: 2,2
#include <m/mat_impl.h>

#undef M
#undef N

#define N 2

// NOTE: 2,M
#include <m/mat_impl.h>

#undef N

#define M 2

// NOTE: N,2
#include <m/mat_impl.h>

#undef M

// NOTE: N,M
#include <m/mat_impl.h>

#undef __m_mat_content_toggle__

#endif
