#ifndef mth_vec_h__
#define mth_vec_h__

/* <mth/vec.h> - vector header
 *      This includes the tvec template class representing a cartesian vector of arbitrary scalar type for dimensions
 *      greater than zero. Included are vector-space operations as well as the inner and outer product. The cross product
 *      and its 2-dimensional equivalent det(a, b) are also defined.
 */

#include <iostream>
#include <iomanip>
#include <array>
#include <type_traits>
#include <cmath>
#include <functional>

#include <mth/mth.h>
#include <mth/comp.h>

namespace mth {

    // Forward declaration to avoid circular includes

    template <typename T, size_t N, size_t M>
    class tmat;

    // N-dimensional vector with scalar type T

    template <typename T, size_t N>
    class tvec {

    private:

        std::array<T, N> values;

    public:

        constexpr tvec() noexcept 
            :values{} {}

        constexpr tvec(const std::array<T, N> &values) noexcept
            :values(values) {}

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N, int>::type = 0>
        constexpr tvec(Q... args) noexcept
            :values{static_cast<T>(args)...} {}

        constexpr operator std::array<T, N>() const noexcept {

            return values;
        }

        template <typename U>
        constexpr operator tvec<U, N>() const noexcept {

            tvec<U, N> result;

            for (size_t i = 0; i < N; i++) {

                result.get(i) = static_cast<U>(get(i));
            }

            return result;
        }

        // Iterator functions

        constexpr auto begin() noexcept {

            return values.begin();
        }

        constexpr auto begin() const noexcept {

            return cbegin();
        }

        constexpr auto end() noexcept {

            return values.end();
        }

        constexpr auto end() const noexcept {

            return cend();
        }

        constexpr auto cbegin() const noexcept {

            return values.cbegin();
        }

        constexpr auto cend() const noexcept {

            return values.cend();
        }

        constexpr auto rbegin() noexcept {

            return values.rbegin();
        }

        constexpr auto rend() noexcept {

            return values.rend();
        }

        constexpr auto crbegin() const noexcept {

            return values.crbegin();
        }

        constexpr auto crend() const noexcept {

            return values.crend();
        }

        constexpr size_t size() const noexcept {

            return N;
        }

        constexpr T &get(size_t index) {

            return values.at(index);
        }

        constexpr const T &get(size_t index) const {

            return values.at(index);
        }

        constexpr T &operator[](size_t index) noexcept {

            return values[index];
        }

        constexpr const T &operator[](size_t index) const noexcept {

            return values[index];
        }

        // Elements for N <= 4 are named

#define BINDING(name, index)    template <size_t n = N, typename std::enable_if<(n == N) && (n > index), int>::type = 0> \
                                constexpr const T & name () const noexcept { return get(index); } \
                                template <size_t n = N, typename std::enable_if<(n == N) && (n > index), int>::type = 0> \
                                constexpr       T & name ()       noexcept { return get(index); }

        BINDING(x, 0)
        BINDING(y, 1)
        BINDING(z, 2)
        BINDING(w, 3)

#undef BINDING

        // TODO: Re-ordered swizzling (maybe find a better way to do this)
        // TODO: Look into returning references here

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 1), int>::type = 0>
        constexpr tvec<T, 2> xy() const noexcept {

            return tvec<T, 2>(get(0), get(1));
        }

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 2), int>::type = 0>
        constexpr tvec<T, 2> yz() const noexcept {

            return tvec<T, 2>(get(1), get(2));
        }

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        constexpr tvec<T, 2> zw() const noexcept {

            return tvec<T, 2>(get(2), get(3));
        }
        
        template <size_t n = N, typename std::enable_if<(n == N) && (n > 2), int>::type = 0>
        constexpr tvec<T, 3> xyz() const noexcept {

            return tvec<T, 3>(get(0), get(1), get(2));
        }

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        constexpr tvec<T, 3> yzw() const noexcept {

            return tvec<T, 3>(get(1), get(2), get(3));
        }

        template <size_t n = N, typename std::enable_if<(n == N) && (n > 3), int>::type = 0>
        constexpr tvec<T, 4> xyzw() const noexcept {

            return tvec<T, 4>(get(0), get(1), get(2), get(3));
        }

        constexpr T magnSqr() const noexcept {

            auto result = static_cast<T>(0);

            for (const auto &value : *this) {

                result += value * value;
            }

            return result;
        }

        template <typename F>
        constexpr auto map(const F &functor) const noexcept {

            using ReturnType = decltype(functor(get(0)));
            tvec<ReturnType, N> result;

            for (size_t i = 0; i < N; i++) {

                result.get(i) = functor(get(i));
            }

            return result;
        }

        // Type is cast to double for more accurate square rooting
        constexpr double magn() const noexcept {

            auto ls = static_cast<double>(magnSqr());

            return std::sqrt(ls);
        }
        
        constexpr T dot(const tvec<T, N> &rhs) const noexcept {

            T result = 0;

            for (size_t i = 0; i < N; i++) {

                result += get(i) * rhs[i];
            }

            return result;
        }

        // Returns a unit vector in the same direction as *this
        constexpr tvec<T, N> unit() const noexcept {

            auto l = static_cast<T>(magn());

            return *this / l;
        }

        constexpr tvec<T, N> &operator+=(const tvec<T, N> &rhs) noexcept {

            for (size_t i = 0; i < N; i++) {

                this->get(i) += rhs[i];
            }

            return *this;
        }

        constexpr tvec<T, N> &operator-=(const tvec<T, N> &rhs) noexcept {

            for (size_t i = 0; i < N; i++) {

                this->get(i) -= rhs[i];
            }

            return *this;
        }

        constexpr tvec<T, N> &operator*=(const T &rhs) noexcept {

            for (size_t i = 0; i < N; i++) {

                this->get(i) *= rhs;
            }

            return *this;
        }

        constexpr tvec<T, N> &operator/=(const T &rhs) noexcept {

            for (size_t i = 0; i < N; i++) {

                this->get(i) /= rhs;
            }

            return *this;
        }
    };

    template <typename T, size_t N>
    constexpr tvec<T, N> operator+(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

        auto result = lhs;

        return result += rhs;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator-(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

        auto result = lhs;

        return result -= rhs;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator-(const tvec<T, N> &rhs) noexcept {

        auto result = rhs;

        for (size_t i = 0; i < N; i++) {

            result[i] = -rhs[i];
        }

        return result;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator*(const T &lhs, const tvec<T, N> &rhs) noexcept {

        auto result = rhs;

        return result *= lhs;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator*(const tvec<T, N> &lhs, const T &rhs) noexcept {

        return rhs * lhs;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator/(const tvec<T, N> &lhs, const T &rhs) noexcept {

        tvec<T, N> result = lhs;

        return result /= rhs;
    }

    template <typename T, size_t N>
    constexpr bool operator==(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

        for (size_t i = 0; i < N; i++) {

            if (!util::isEqual(lhs[i], rhs[i])) return false;
        }

        return true;
    }

    template <typename T, size_t N>
    constexpr bool operator!=(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T, size_t N>
    std::ostream &operator<<(std::ostream &lhs, const tvec<T, N> &rhs) {

        lhs << "(";

        for (size_t i = 0; i < N - 1; i++) {

            lhs << rhs[i] << ", ";
        }

        return lhs << rhs[N - 1] << ")";
    }

    // Specialized "static" functions

    namespace vec {

        template <typename T>
        constexpr tvec<T, 3> cross(const tvec<T, 3> &lhs, const tvec<T, 3> &rhs) noexcept {

            return tvec<T, 3>(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                              lhs[2] * rhs[0] - lhs[0] * rhs[2],
                              lhs[0] * rhs[1] - lhs[1] * rhs[0]);
        }

        template <typename T>
        T constexpr det(const tvec<T, 2> &lhs, const tvec<T, 2> &rhs) noexcept {

            return lhs[0] * rhs[1] - lhs[1] * rhs[0]; // tmat<T, 2, 2>(lhs, rhs).det()
        }
    }

    // Alias types for single-digit dimension and scalar type int, long, float, double and mth::comp

#define CREATE_ALIASES(n) using ivec ## n = tvec<int, n>; \
                          using lvec ## n = tvec<long, n>; \
                          using fvec ## n = tvec<float, n>; \
                          using dvec ## n = tvec<double, n>; \
                          using cvec ## n = tvec<mth::comp, n>; \
                          using vec ## n = dvec ## n;

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

    namespace vec {

        /**
         *
         * Apply a functor element-wise over multiple vectors, the multiple parameter version of
         * tvec::map.
         *
         * `map(functor, vec0, vec1, ...)`
         *
         * is loosely equivalent to:
         *
         * `tvec{functor(vec0.get(0), vec1.get(0), ...),
         *       functor(vec0.get(1), vec1.get(1), ...),
         *       ...}`
         *
         */
        template <typename F, typename... Args, size_t N>
        constexpr auto map(F functor, tvec<Args, N>... args) {

            using ReturnType = decltype(functor(args.get(0)...));
            tvec<ReturnType, N> result;

            for (size_t i = 0; i < N; i++) {

                result.get(i) = functor(args.get(i)...);
            }

            return result;
        }

        template <typename T, size_t N>
        constexpr tvec<T, N> hadamard(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

            tvec<T, N> result;

            for (size_t i = 0; i < N; i++) {

                result.get(i) = lhs.get(i) * rhs.get(i);
            }

            return result;
        }

        template <typename T, size_t N>
        constexpr tmat<T, N, N> outerProduct(const tvec<T, N> &lhs, const tvec<T, N> &rhs) noexcept {

            auto lMat = mth::tmat<T, 1, N>(lhs);
            auto rMatT = mth::tmat<T, N, 1>(rhs);

            return lMat * rMatT;
        }

        template <typename T, size_t N>
        constexpr T dot(const mth::tvec<T, N> &lhs, const mth::tvec<T, N> &rhs) {

            return lhs.dot(rhs);
        }
    }

    // Constants for the axes in 3-dimensions

    template <typename T>
    constexpr tvec<T, 3> X_AXIS = tvec<T, 3>(1, 0, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Y_AXIS = tvec<T, 3>(0, 1, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Z_AXIS = tvec<T, 3>(0, 0, 1);
    
    template <typename T, size_t N>
    constexpr double abs(const mth::tvec<T, N> &x) noexcept {

        return x.magn();
    }
}

namespace std {

    // Hash operator for use in certain STL containers
    template <typename T, size_t N>
    struct hash<mth::tvec<T, N>> {

        size_t operator()(const mth::tvec<T, N> &x) {

            size_t result = 0;

            for (size_t i = 0; i < N; i++) {

                result ^= hash<T>()(x[i]);
            }

            return result;
        }
    };
}

#endif
