#ifndef __Maths_tmat2_h__
#define __Maths_tmat2_h__

#include <limits>
#include <cmath>
#include <ostream>

#ifndef __Maths_vec2_h__
#include <Maths/vec2.h>
#endif

#define N 2

namespace m {

    // requires 0, 1, T + T, T - T, T * T, T / T, std::abs(T), std::numeric_limits<T>::epsilon(), std::ostream << T
    template <typename T> 
    struct tmat2 {

        T values[N * N]; // column-major

        static size_t getIndex(size_t x, size_t y);

        const T &get(size_t x, size_t y) const;
        T &get(size_t x, size_t y);
        tvec2<T> getRow(size_t y) const;
        tvec2<T> getColumn(size_t x) const;

        tmat2() :tmat2{1, 0, 0, 1} {}
        tmat2(tvec2<T> c1, tvec2<T> c2) :tmat2{c1.x, c1.y, c2.x, c2.y} {}
        tmat2(T values[N * N]) :tmat2{values[0], values[1], values[2], values[3]} {}
        tmat2(std::initializer_list<T> values);

        T det() const;
        tmat2 cofactors() const;
        tmat2 transpose() const;
        tmat2 adjoint() const;
        tmat2 inverse() const;
        tmat2 unit() const;

        inline static tmat2 identity() { return tmat2(); }
    };

    typedef tmat2<int>    imat2;
    typedef tmat2<long>   lmat2;
    typedef tmat2<float>   mat2;
    typedef tmat2<double> dmat2;

    template <typename T>
    tmat2<T> operator+(const tmat2<T> &a, const tmat2<T> &b);

    template <typename T>
    tmat2<T> operator-(const tmat2<T> &a, const tmat2<T> &b);

    template <typename T>
    tmat2<T> operator*(T scalar, const tmat2<T> &matrix);

    template <typename T>
    tmat2<T> operator*(const tmat2<T> &matrix, T scalar);

    template <typename T>
    tmat2<T> operator/(const tmat2<T> &matrix, T scalar);

    template <typename T>
    tmat2<T> operator*(const tmat2<T> &a, const tmat2<T> &b);

    template <typename T>
    tmat2<T> operator/(const tmat2<T> &a, const tmat2<T> &b);

    template <typename T>
    tvec2<T> operator*(const tmat2<T> &matrix, const tvec2<T> &vector);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tmat2<T> &matrix);
}

// Template implementation

template <typename T>
m::tmat2<T>::tmat2(std::initializer_list<T> init) {

    size_t i = 0;

    for (auto v : init) {

        if (i > N * N - 1) break;
        values[i++] = v;
    }
}

template <typename T>
const T &m::tmat2<T>::get(size_t x, size_t y) const {

    return values[getIndex(x, y)];
}

template <typename T>
T &m::tmat2<T>::get(size_t x, size_t y) {

    return values[getIndex(x, y)];
}

template <typename T>
size_t m::tmat2<T>::getIndex(size_t x, size_t y) {

    if (x > N - 1) throw std::invalid_argument("mat2 only has 2 columns!");
    if (y > N - 1) throw std::invalid_argument("mat2 only has 2 rows!");

    return x * N + y;
}

template <typename T>
m::tvec2<T> m::tmat2<T>::getRow(size_t y) const {

    return tvec2<T>(get(0, y), get(1, y));
}

template <typename T>
m::tvec2<T> m::tmat2<T>::getColumn(size_t x) const {

    return tvec2<T>(get(x, 0), get(x, 1));
}

template <typename T>
T m::tmat2<T>::det() const {

    return values[0] * values[3] - values[2] * values[1];
}

template <typename T>
m::tmat2<T> m::tmat2<T>::cofactors() const {

    return tmat2<T>{values[3], -values[2], -values[1], values[0]};
}

template <typename T>
m::tmat2<T> m::tmat2<T>::transpose() const {

    return tmat2<T>{values[0], values[2], values[1], values[3]};
}

template <typename T>
m::tmat2<T> m::tmat2<T>::adjoint() const {

    return cofactors().transpose();
}

template <typename T>
m::tmat2<T> m::tmat2<T>::inverse() const {

    T determinant = det();

    if (util::checkZero(determinant)) throw std::invalid_argument("Singular matrix can't be inverted");

    tmat2<T> adj = adjoint();

    return adj / determinant;
}

template <typename T>
m::tmat2<T> m::tmat2<T>::unit() const {

    return *this / det();
}

template <typename T>
m::tmat2<T> m::operator+(const m::tmat2<T> &a, const m::tmat2<T> &b) {

    tmat2<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] += b.values[i];

    return result;
}

template <typename T>
m::tmat2<T> m::operator-(const tmat2<T> &a, const m::tmat2<T> &b) {

    tmat2<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] -= b.values[i];

    return result;
}

template <typename T>
m::tmat2<T> m::operator*(T scalar, const m::tmat2<T> &matrix) {

    tmat2<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] *= scalar;

    return result;
}

template <typename T>
m::tmat2<T> m::operator*(const m::tmat2<T> &matrix, T scalar) {

    return scalar * matrix;
}

template <typename T>
m::tmat2<T> m::operator/(const m::tmat2<T> &matrix, T scalar) {

    tmat2<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] /= scalar;

    return result;
}

template <typename T>
m::tmat2<T> m::operator*(const m::tmat2<T> &a, const tmat2<T> &b) {

    tmat2<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = tvec2<T>::dot(a.getRow(y), b.getColumn(x));
        }
    }

    return result;
}

template <typename T>
m::tmat2<T> m::operator/(const m::tmat2<T> &a, const m::tmat2<T> &b) {

    return b.inverse() * a;
}

template <typename T>
m::tvec2<T> m::operator*(const m::tmat2<T> &matrix, const m::tvec2<T> &vector) {

    tvec2<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(y) = tvec2<T>::dot(matrix.getRow(y), vector);
        }
    }

    return result;
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const m::tmat2<T> &matrix) {

    return stream << std::endl << "__\t\t\t\b \b__" << std::endl
                               << "|\t" << matrix.values[0] << "\t" << matrix.values[2] << "\t|" << std::endl
                               << "|\t" << matrix.values[1] << "\t" << matrix.values[3] << "\t|" << std::endl
                               << "--\t\t\t\b \b--";
}

#undef N

#endif
