#ifndef __m_vec_h__
#define __m_vec_h__

#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <array>
#include <type_traits>

#ifndef __m_constants_h__
#include <m/constants.h>
#endif

namespace m {

    template <typename T, size_t N>
    class tvec {

    private:

        std::array<T, N> values;

    public:

        constexpr tvec() noexcept {

            values.fill(0);
        }

        tvec(const std::array<T, N> &values) noexcept
            :values(values) {}

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N, int>::type = 0>
        constexpr tvec(Q... args) noexcept
            :values{static_cast<T>(args)...} {}

        operator std::array<T, N>() noexcept {

            return values;
        }

        auto begin() {

            return values.begin();
        }

        auto begin() const {

            return cbegin();
        }

        auto end() {

            return values.end();
        }

        auto end() const {

            return cend();
        }

        auto cbegin() const {

            return values.cbegin();
        }

        auto cend() const {

            return values.cend();
        }

        auto rbegin() { 

            return values.rbegin();
        }

        auto rend() {

            return values.rend();
        }

        auto crbegin() const {

            return values.crbegin();
        }

        auto crend() const {

            return values.crend();
        }

        constexpr size_t size() const noexcept {

            return N;
        }

        auto &get(size_t index) {

            if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

            return values[index];
        }

        const auto &get(size_t index) const {

            if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

            return values[index];
        }

#define BINDING(name, value) const auto & name () const { return value ; } \
                                   auto & name ()       { return value; }

        BINDING(x, this->get(0))
        BINDING(y, this->get(1))
        BINDING(z, this->get(2))
        BINDING(w, this->get(3))

#undef BINDING

        // TODO: There should be a better way of doing this

        tvec<T, 2> xy() const {

            return tvec<T, 2>(x(), y());
        }

        tvec<T, 2> yz() const {

            return tvec<T, 2>(y(), z());
        }

        tvec<T, 2> zw() const {

            return tvec<T, 2>(z(), w());
        }
        
        tvec<T, 3> xyz() const {

            return tvec<T, 3>(x(), y(), z());
        }

        tvec<T, 3> yzw() const {

            return tvec<T, 3>(y(), z(), w());
        }

        auto magnSqr() const noexcept {

            auto result = static_cast<T>(0);

            for (const auto &value : *this) {

                result += value * value;
            }

            return result;
        }

        auto magn() const noexcept {

            auto ls = static_cast<double>(magnSqr());

            if (util::checkZero(ls)) return 0.0;

            return std::sqrt(ls);
        }
        
        auto dot(const tvec<T, N> &rhs) const {

            T result = 0;

            for (size_t i = 0; i < N; i++) {

                result += this->get(i) * rhs.get(i);
            }

            return result;
        }

        auto unit() const {

            auto l = static_cast<T>(magn());

            if (util::checkZero(l)) throw std::invalid_argument("m::exception: unit() called on zero vector");

            return *this / l;
        }

        auto &operator+=(const tvec<T, N> &rhs) {

            for (size_t i = 0; i < N; i++) {

                this->get(i) += rhs.get(i);
            }

            return *this;
        }

        auto &operator-=(const tvec<T, N> &rhs) {

            for (size_t i = 0; i < N; i++) {

                this->get(i) -= rhs.get(i);
            }

            return *this;
        }

        auto &operator*=(const T &rhs) {

            for (size_t i = 0; i < N; i++) {

                this->get(i) *= rhs;
            }

            return *this;
        }

        auto &operator/=(const T &rhs) {

            for (size_t i = 0; i < N; i++) {

                this->get(i) /= rhs;
            }

            return *this;
        }

        friend auto operator+(const tvec<T, N> &lhs, const tvec<T, N> &rhs) {

            auto result = lhs;

            return result += rhs;
        }

        friend auto operator-(const tvec<T, N> &lhs, const tvec<T, N> &rhs) {

            auto result = lhs;

            return result -= rhs;
        }

        friend auto operator-(const tvec<T, N> &rhs) {

            auto result = rhs;

            for (size_t i = 0; i < N; i++) {

                result.get(i) = -rhs.get(i);
            }

            return result;
        }

        friend auto operator*(const T &lhs, const tvec<T, N> &rhs) {

            auto result = rhs;

            return result *= lhs;
        }

        friend auto operator*(const tvec<T, N> &lhs, const T &rhs) {

            return rhs * lhs;
        }

        friend auto operator/(const tvec<T, N> &lhs, const T &rhs) {

            tvec<T, N> result = lhs;

            return result /= rhs;
        }

        friend auto &operator<<(std::ostream &lhs, const tvec<T, N> &rhs) {

            lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

            for (size_t i = 0; i < N - 1; i++) {

                lhs << rhs.get(i) << ", ";
            }

            return lhs << rhs.get(N - 1) << ")";
        }
    };

    // NOTE: Static functions for vectors

    namespace vec {

        template <typename T, size_t N>
        auto dot(const tvec<T, N> &lhs, const tvec<T, N> &rhs) {

            return lhs.dot(rhs);
        }

        template <typename T>
        auto cross(const tvec<T, 3> &lhs, const tvec<T, 3> &rhs) {

            return tvec<T, 3>(lhs.get(2) * rhs.get(3) - lhs.get(3) * rhs.get(2),
                              lhs.get(3) * rhs.get(1) - lhs.get(1) * rhs.get(3),
                              lhs.get(1) * rhs.get(2) - lhs.get(2) * rhs.get(1));
        }
    }

    // TODO: Maybe move to double as default
    
#define CREATE_ALIASES(n) using ivec ## n = tvec<int, n>; \
                          using lvec ## n = tvec<long, n>; \
                          using  vec ## n = tvec<float, n>; \
                          using dvec ## n = tvec<double, n>; \
                          using cvec ## n = tvec<std::complex<double>, n>;

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

#endif
