#ifndef __m_vec_h__
#define __m_vec_h__

/* <m/vec.h> - vector header
 *      This includes the tvec template class representing a cartesian vector of arbitrary scalar type for dimensions
 *      greater than zero. Included are vector-space operations as well as the inner and outer product. The cross product
 *      and its 2-dimensional equivalent det(a, b) are also defined.
 */

#include <iostream>
#include <iomanip>
#include <array>
#include <type_traits>

#include <m/comp.h>

namespace m {

    // Forward declaration to avoid circular includes

    template <typename T, size_t N, size_t M>
    class tmat;

    // Forward declaration for friending

    template <typename T, size_t N>
    class tvec;

    template <typename T, size_t N>
    tvec<T, N> operator+(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    tvec<T, N> operator-(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    tvec<T, N> operator-(const tvec<T, N> &rhs);

    template <typename T, size_t N>
    tvec<T, N> operator*(const T &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    tvec<T, N> operator*(const tvec<T, N> &lhs, const T &rhs);

    template <typename T, size_t N>
    tvec<T, N> operator/(const tvec<T, N> &lhs, const T &rhs);

    template <typename T, size_t N>
    bool operator==(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    bool operator!=(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    std::ostream &operator<<(std::ostream &lhs, const tvec<T, N> &rhs);

    // TODO: noexcept on operators?

    // N-dimensional vector with scalar type T

    template <typename T, size_t N>
    class tvec {

    private:

        std::array<T, N> values;

    public:

        constexpr tvec() noexcept {

            values.fill(0);
        }

        tvec(const std::array<T, N> &values) noexcept;

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N, int>::type = 0>
        constexpr tvec(Q... args) noexcept
            :values{static_cast<T>(args)...} {}

        constexpr operator std::array<T, N>() noexcept {

            return values;
        }

        // Iterator functions

        auto begin();
        auto begin() const;
        auto end();
        auto end() const;

        auto cbegin() const;
        auto cend() const;

        auto rbegin();
        auto rend();

        auto crbegin() const;
        auto crend() const;

        constexpr size_t size() const noexcept {

            return N;
        }

        T &get(size_t index);
        const T &get(size_t index) const;

        // Elements for N <= 4 are named

#define BINDING(name, index)    template <size_t n = N, typename std::enable_if<(n == N) && (n > index), int>::type = 0> \
                                const T & name () const; \
                                template <size_t n = N, typename std::enable_if<(n == N) && (n > index), int>::type = 0> \
                                      T & name ();

        BINDING(x, 0)
        BINDING(y, 1)
        BINDING(z, 2)
        BINDING(w, 3)

#undef BINDING

        // TODO: Re-ordered swizzling (maybe find a better way to do this)
        // TODO: Look into returning references here

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 1), int>::type = 0>
        tvec<T, 2> xy() const;

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 2), int>::type = 0>
        tvec<T, 2> yz() const;

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        tvec<T, 2> zw() const;
        
        template <size_t n = N, typename std::enable_if<(n == N) && (n > 2), int>::type = 0>
        tvec<T, 3> xyz() const;

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        tvec<T, 3> yzw() const;

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        tvec<T, 4> xyzw() const;

        T magnSqr() const noexcept;

        // Type is cast to double for more accurate square rooting
        double magn() const noexcept;
        
        T dot(const tvec<T, N> &rhs) const;
        static T dot(const tvec<T, N> &lhs, const tvec<T, N> &rhs);
        static tmat<T, N, N> outerProduct(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

        // Returns a unit vector in the same direction as *this
        tvec<T, N> unit() const;

        tvec<T, N> &operator+=(const tvec<T, N> &rhs);
        tvec<T, N> &operator-=(const tvec<T, N> &rhs);
        tvec<T, N> &operator*=(const T &rhs);
        tvec<T, N> &operator/=(const T &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator+(const tvec<U, M> &lhs, const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator-(const tvec<U, M> &lhs, const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator-(const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator*(const U &lhs, const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator*(const tvec<U, M> &lhs, const U &rhs);

        template <typename U, size_t M>
        friend tvec<U, M> operator/(const tvec<U, M> &lhs, const U &rhs);

        template <typename U, size_t M>
        friend bool operator==(const tvec<U, M> &lhs, const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend bool operator!=(const tvec<U, M> &lhs, const tvec<U, M> &rhs);

        template <typename U, size_t M>
        friend std::ostream &operator<<(std::ostream &lhs, const tvec<U, M> &rhs);
    };

    // Specialized "static" functions

    namespace vec {

        template <typename T>
        tvec<T, 3> cross(const tvec<T, 3> &lhs, const tvec<T, 3> &rhs);

        template <typename T>
        T det(const tvec<T, 2> &lhs, const tvec<T, 2> &rhs);
    }

    // Alias types for single-digit dimension and scalar type int, long, float, double and m::comp

    // TODO: Maybe move to double as default
#define CREATE_ALIASES(n) using ivec ## n = tvec<int, n>; \
                          using lvec ## n = tvec<long, n>; \
                          using  vec ## n = tvec<float, n>; \
                          using dvec ## n = tvec<double, n>; \
                          using cvec ## n = tvec<m::comp, n>;

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

    // Constants for the axes in 3-dimensions

    template <typename T>
    constexpr tvec<T, 3> X_AXIS = tvec<T, 3>(1, 0, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Y_AXIS = tvec<T, 3>(0, 1, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Z_AXIS = tvec<T, 3>(0, 0, 1);
}

#include <m/vec_impl.h>

#endif
