#ifndef __Maths_mat3_h__
#define __Maths_mat3_h__

#include <limits>
#include <cmath>
#include <ostream>

#ifndef __Maths_vec3_h__
#include <Maths/vec3.h>
#endif

#ifndef __Maths_mat2_h__
#include <Maths/mat2.h>
#endif

#define N 3

namespace m {

    template <typename T> // requires 0, 1, -1, T + T, T - T, T * T, T / T, std::abs(T), std::numeric_limits<T>::epsilon(), std::ostream << T
    struct tmat3 {

        T values[N * N]; // column-major

        static size_t getIndex(size_t x, size_t y);

        const T &get(size_t x, size_t y) const;
        T &get(size_t x, size_t y);
        tvec3<T> getRow(size_t y) const;
        tvec3<T> getColumn(size_t x) const;

        tmat3() :tmat3{1, 0, 0, 0, 1, 0, 0, 0, 1} {}
        tmat3(tvec3<T> c1, tvec3<T> c2, tvec3<T> c3) :tmat3{c1.x, c1.y, c1.z, c2.x, c2.y, c2.z, c3.x, c3.y, c3.z} {}
        tmat3(T values[N * N]) :tmat3{values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]} {}
        tmat3(std::initializer_list<T> values);

        T det() const;
        tmat2<T> getMinor(size_t x, size_t y) const;
        tmat3 cofactors() const;
        tmat3 transpose() const;
        tmat3 adjoint() const;
        tmat3 inverse() const;

        static tmat3 identity();
        static tmat3 scale(T x, T y, T z);
    };

    typedef tmat3<int>    imat3;
    typedef tmat3<long>   lmat3;
    typedef tmat3<float>   mat3;
    typedef tmat3<double> dmat3;

    template <typename T>
    tmat3<T> operator+(const tmat3<T> &a, const tmat3<T> &b);

    template <typename T>
    tmat3<T> operator-(const tmat3<T> &a, const tmat3<T> &b);

    template <typename T>
    tmat3<T> operator*(T scalar, const tmat3<T> &matrix);

    template <typename T>
    tmat3<T> operator*(const tmat3<T> &matrix, T scalar);

    template <typename T>
    tmat3<T> operator/(const tmat3<T> &matrix, T scalar);

    template <typename T>
    tmat3<T> operator*(const tmat3<T> &a, const tmat3<T> &b);

    template <typename T>
    tmat3<T> operator/(const tmat3<T> &a, const tmat3<T> &b);

    template <typename T>
    tvec3<T> operator*(const tmat3<T> &matrix, const tvec3<T> &vector);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tmat3<T> &matrix);
}

// Template implementation

template <typename T>
m::tmat3<T>::tmat3(std::initializer_list<T> init) {

    size_t i = 0;

    for (auto v : init) {

        if (i > N * N - 1) break;
        values[i++] = v;
    }
}

template <typename T>
size_t m::tmat3<T>::getIndex(size_t x, size_t y) {

    if (x > N - 1) throw std::invalid_argument("mat3 only has 3 columns!");
    if (y > N - 1) throw std::invalid_argument("mat3 only has 3 rows!");

    return x * N + y;
}

template <typename T>
const T &m::tmat3<T>::get(size_t x, size_t y) const {

    return values[getIndex(x, y)];
}

template <typename T>
T &m::tmat3<T>::get(size_t x, size_t y) {

    return values[getIndex(x, y)];
}

template <typename T>
m::tvec3<T> m::tmat3<T>::getRow(size_t y) const {

    return tvec3<T>(get(0, y), get(1, y), get(2, y));
}

template <typename T>
m::tvec3<T> m::tmat3<T>::getColumn(size_t x) const {

    return tvec3<T>(get(x, 0), get(x, 1), get(x, 2));
}

template <typename T>
m::tmat2<T> m::tmat3<T>::getMinor(size_t x, size_t y) const {

    size_t rx = 0;
    size_t ry = 0;

    tmat2<T> result;

    for (size_t ix = 0; ix < N; ix++) {

        if (ix == x) continue;

        for (size_t iy = 0; iy < N; iy++) {

            if (iy == y) continue;

            result.get(rx, ry) = get(ix, iy);

            ry++;
        }

        rx++;
        ry = 0;
    }

    return result;
}

template <typename T>
T m::tmat3<T>::det() const {

    T result = 0;
    T s = 1;
    size_t row = 0; // can be any row

    for (size_t x = 0; x < N; x++) {

        result += s * get(x, row) * getMinor(x, row).det();
        s *= -1;
    }

    return result;
}

template <typename T>
m::tmat3<T> m::tmat3<T>::cofactors() const {

    tmat3<T> result;
    T s = 1;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = s * getMinor(x, y).det();
            s *= -1;
        }
    }

    return result;
}

template <typename T>
m::tmat3<T> m::tmat3<T>::transpose() const {

    tmat3<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = get(y, x);
        }
    }

    return result;
}

template <typename T>
m::tmat3<T> m::tmat3<T>::adjoint() const {

    return cofactors().transpose();
}

template <typename T>
m::tmat3<T> m::tmat3<T>::inverse() const {

    T determinant = det();

    if (util::checkZero(determinant)) throw std::invalid_argument("Singular matrix can't be inverted");

    tmat3<T> adj = adjoint();

    return adj / determinant;
}

template <typename T>
m::tmat3<T> m::tmat3<T>::identity() {

    return tmat3<T>();
}

template <typename T>
m::tmat3<T> m::tmat3<T>::scale(T x, T y, T z) {

    return tmat3<T>{x, 0, 0, 0, y, 0, 0, 0, z};
}

template <typename T>
m::tmat3<T> m::operator+(const m::tmat3<T> &a, const m::tmat3<T> &b) {

    tmat3<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] += b.values[i];

    return result;
}

template <typename T>
m::tmat3<T> m::operator-(const tmat3<T> &a, const m::tmat3<T> &b) {

    tmat3<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] -= b.values[i];

    return result;
}

template <typename T>
m::tmat3<T> m::operator*(T scalar, const m::tmat3<T> &matrix) {

    tmat3<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] *= scalar;

    return result;
}

template <typename T>
m::tmat3<T> m::operator*(const m::tmat3<T> &matrix, T scalar) {

    return scalar * matrix;
}

template <typename T>
m::tmat3<T> m::operator/(const m::tmat3<T> &matrix, T scalar) {

    tmat3<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] /= scalar;

    return result;
}

template <typename T>
m::tmat3<T> m::operator*(const m::tmat3<T> &a, const tmat3<T> &b) {

    tmat3<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = tvec3<T>::dot(a.getRow(y), b.getColumn(x));
        }
    }

    return result;
}

template <typename T>
m::tmat3<T> m::operator/(const m::tmat3<T> &a, const m::tmat3<T> &b) {

    return b.inverse() * a;
}

template <typename T>
m::tvec3<T> m::operator*(const m::tmat3<T> &matrix, const m::tvec3<T> &vector) {

    tvec3<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(y) = tvec3<T>::dot(matrix.getRow(y), vector);
        }
    }

    return result;
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const m::tmat3<T> &matrix) {

    return stream << std::endl << "__\t\t\t\t\b \b__" << std::endl
                               << "|\t" << matrix.values[0] << "\t" << matrix.values[3] << "\t" << matrix.values[6]<< "\t|" << std::endl
                               << "|\t" << matrix.values[1] << "\t" << matrix.values[4] << "\t" << matrix.values[7]<< "\t|" << std::endl
                               << "|\t" << matrix.values[2] << "\t" << matrix.values[5] << "\t" << matrix.values[8]<< "\t|" << std::endl
                               << "--\t\t\t\t\b \b--";
}

#undef N

#endif
