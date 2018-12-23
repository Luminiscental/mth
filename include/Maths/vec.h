#ifndef __Maths_vec_h__
#define __Maths_vec_h__

#include <memory>
#include <complex>
#include <cmath>
#include <ostream>

#ifndef __Maths_constants_h__
#include <Maths/constants.h>
#endif

namespace m {

    template <typename T, size_t N>
    struct tvec {

        T values[N * N];

        tvec() {

            for (size_t i = 0; i < N; i++) {

                values[i] = 0;
            }
        }

        tvec(const tvec<T, N> &other) {

            for (size_t i = 0; i < N; i++) {

                values[i] = other.get(i);
            }
        }

        tvec(std::initializer_list<T> init) {
                                                              
            size_t i = 0;
                                                              
            for (T value : init) {
                                                              
                values[i++] = value;
            }
        }

        T &get(size_t index) {

            return values[index];
        }

        const T &get(size_t index) const {

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

            return std::sqrt(static_cast<double>(magnSqr()));
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

            if (util::checkZero(l)) throw std::invalid_argument("zero vector has no unit equivalent");

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

            stream << "(";

            for (size_t i = 0; i < N - 1; i++) {

                stream << vector.get(i) << ", ";
            }

            return stream << vector.get(N - 1) << ")";
        }
    };

    template <typename T>
    tvec<T, 3> x_axis{1, 0, 0};

    template <typename T>
    tvec<T, 3> y_axis{0, 1, 0};

    template <typename T>
    tvec<T, 3> z_axis{0, 0, 1};
    
#define TYPEDEF_VEC(n) typedef tvec<int, n>    ivec ## n; \
                       typedef tvec<long, n>   lvec ## n; \
                       typedef tvec<float, n>   vec ## n; \
                       typedef tvec<double, n> dvec ## n; \
                       typedef tvec<std::complex<double>, n> cvec ## n;

    TYPEDEF_VEC(2)
    TYPEDEF_VEC(3)
    TYPEDEF_VEC(4)

#undef TYPEDEF_VEC

}

#endif
