#ifndef __Maths_vec_h__
#define __Maths_vec_h__

#include <complex>
#include <cmath>

#include <iostream>
#include <iomanip>

#include <array>
#include <type_traits>

#ifndef __Maths_constants_h__
#include <Maths/constants.h>
#endif

namespace m {

    template <typename T, size_t N>
    class tvec {

    private:
        std::array<T, N> values;

    public:
        tvec() {

            values.fill(0);
        }

        tvec(const std::array<T, N> &values) : values(values) {} 

        tvec(const tvec<T, N> &other) : tvec(other.values) {}

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N, int>::type = 0>
        tvec(Q... args) : values{static_cast<T>(args)...} {}

        T &get(size_t index) {

            if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

            return values[index];
        }

        const T &get(size_t index) const {

            if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

            return values[index];
        }

        T x() const {

            return get(0);
        }

        T y() const {

            return get(1);
        }

        T z() const {

            return get(2);
        }
        
        T w() const {

            return get(3);
        }

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

        T magnSqr() const {

            T result = 0;

            for (size_t i = 0; i < N; i++) {

                result += this->get(i) * this->get(i);
            }

            return result;
        }

        double magn() const {

            double ls = static_cast<double>(magnSqr());

            if (util::checkZero(ls)) return 0.0;

            return std::sqrt(ls);
        }
        
        T dot(const tvec<T, N> &other) {

            T result = 0;

            for (size_t i = 0; i < N; i++) {

                result += this->get(i) * other.get(i);
            }

            return result;
        }

        tvec<T, N> unit() const {

            T l = static_cast<T>(magn());

            if (util::checkZero(l)) throw std::invalid_argument("m::exception: unit() called on zero vector");

            return *this / l;
        }

        friend tvec<T, N> operator+(const tvec<T, N> &a, const tvec<T, N> &b) {

            tvec<T, N> result = a;

            for (size_t i = 0; i < N; i++) {

                result.get(i) += b.get(i);
            }

            return result;
        }

        friend tvec<T, N> operator-(const tvec<T, N> &a, const tvec<T, N> &b) {

            tvec<T, N> result = a;

            for (size_t i = 0; i < N; i++) {

                result.get(i) -= b.get(i);
            }

            return result;
        }

        friend tvec<T, N> operator-(const tvec<T, N> &vector) {

            tvec<T, N> result = vector;

            for (size_t i = 0; i < N; i++) {

                vector.get(i) = -vector.get(i);
            }

            return result;
        }

        friend tvec<T, N> operator*(T scalar, const tvec<T, N> &vector) {

            tvec<T, N> result = vector;

            for (size_t i = 0; i < N; i++) {

                result.get(i) *= scalar;
            }

            return result;
        }

        friend tvec<T, N> operator*(const tvec<T, N> &vector, T scalar) {

            return scalar * vector;
        }

        friend tvec<T, N> operator/(const tvec<T, N> &vector, T scalar) {

            tvec<T, N> result = vector;

            for (size_t i = 0; i < N; i++) {

                result.get(i) /= scalar;
            }

            return result;
        }

        friend std::ostream &operator<<(std::ostream &stream, const tvec<T, N> &vector) {

            stream << std::fixed << std::setprecision(2) << "(";

            for (size_t i = 0; i < N - 1; i++) {

                stream << vector.get(i) << ", ";
            }

            return stream << vector.get(N - 1) << ")";
        }
    };

#define TYPEDEF_VEC(n) typedef tvec<int, n>    ivec ## n; \
                       typedef tvec<long, n>   lvec ## n; \
                       typedef tvec<float, n>   vec ## n; \
                       typedef tvec<double, n> dvec ## n; \
                       typedef tvec<std::complex<double>, n> cvec ## n;

    // TODO: Maybe move to double as default
    
    TYPEDEF_VEC(2)
    TYPEDEF_VEC(3)
    TYPEDEF_VEC(4)

#undef TYPEDEF_VEC

}

#endif
