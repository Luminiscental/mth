#ifndef __Maths_mat4_h__
#define __Maths_mat4_h__

#include <limits>
#include <cmath>
#include <ostream>

#ifndef __Maths_vec4_h__
#include <Maths/vec4.h>
#endif

#ifndef __Maths_vec4_h__
#include <Maths/vec4.h>
#endif

#define N 4

namespace m {

    template <typename T> // requires 0, 1, -1, 2, -2, T + T, T - T, T * T, T / T, static_cast<double>(T), static_cast<T>(double), std::abs(T), std::numeric_limits<T>::epsilon(), std::ostream << T
    struct tmat4 {

        T values[N * N]; // column-major

        static size_t getIndex(size_t x, size_t y);

        const T &get(size_t x, size_t y) const;
        T &get(size_t x, size_t y);
        tvec4<T> getRow(size_t y) const;
        tvec4<T> getColumn(size_t x) const;

        tmat4() :tmat4{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1} {}
        tmat4(tvec4<T> c1, tvec4<T> c2, tvec4<T> c3, tvec4<T> c4) :tmat4{c1.x, c1.y, c1.z, c1.w, c2.x, c2.y, c2.z, c2.w, c3.x, c3.y, c3.z, c3.w, c4.x, c4.y, c4.z, c4.w} {}
        tmat4(T values[N * N]) :tmat4{values[0], values[1], values[2], values[4], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], values[15]} {}
        tmat4(std::initializer_list<T> values);

        T det() const;
        tmat3<T> getMinor(size_t x, size_t y) const;
        tmat4 cofactors() const;
        tmat4 transpose() const;
        tmat4 adjoint() const;
        tmat4 inverse() const;

        static tmat4 identity();
        static tmat4 scale(T x, T y, T z);
        static tmat4 scale(tvec3<T> v);
        static tmat4 translate(T x, T y, T z);
        static tmat4 orthogonalProjection(T left, T right, T bottom, T top, T near, T far);
        static tmat4 perspectiveProjection(T fovHorizontal, T fovVertical, T near, T far);
    };

    typedef tmat4<int>    imat4;
    typedef tmat4<long>   lmat4;
    typedef tmat4<float>   mat4;
    typedef tmat4<double> dmat4;

    template <typename T>
    tmat4<T> operator+(const tmat4<T> &a, const tmat4<T> &b);

    template <typename T>
    tmat4<T> operator-(const tmat4<T> &a, const tmat4<T> &b);

    template <typename T>
    tmat4<T> operator*(T scalar, const tmat4<T> &matrix);

    template <typename T>
    tmat4<T> operator*(const tmat4<T> &matrix, T scalar);

    template <typename T>
    tmat4<T> operator/(const tmat4<T> &matrix, T scalar);

    template <typename T>
    tmat4<T> operator*(const tmat4<T> &a, const tmat4<T> &b);

    template <typename T>
    tmat4<T> operator/(const tmat4<T> &a, const tmat4<T> &b);

    template <typename T>
    tvec4<T> operator*(const tmat4<T> &matrix, const tvec4<T> &vector);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tmat4<T> &matrix);
}

// Template implementation

template <typename T>
m::tmat4<T>::tmat4(std::initializer_list<T> init) {

    size_t i = 0;

    for (auto v : init) {

        if (i > N * N - 1) break;
        values[i++] = v;
    }
}

template <typename T>
size_t m::tmat4<T>::getIndex(size_t x, size_t y) {

    if (x > N - 1) throw std::invalid_argument("mat4 only has 4 columns!");
    if (y > N - 1) throw std::invalid_argument("mat4 only has 4 rows!");

    return x * N + y;
}

template <typename T>
const T &m::tmat4<T>::get(size_t x, size_t y) const {

    return values[getIndex(x, y)];
}

template <typename T>
T &m::tmat4<T>::get(size_t x, size_t y) {

    return values[getIndex(x, y)];
}

template <typename T>
m::tvec4<T> m::tmat4<T>::getRow(size_t y) const {

    return tvec4<T>(get(0, y), get(1, y), get(2, y), get(3, y));
}

template <typename T>
m::tvec4<T> m::tmat4<T>::getColumn(size_t x) const {

    return tvec4<T>(get(x, 0), get(x, 1), get(x, 2), get(x, 3));
}

template <typename T>
m::tmat3<T> m::tmat4<T>::getMinor(size_t x, size_t y) const {

    size_t rx = 0;
    size_t ry = 0;

    tmat3<T> result;

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
T m::tmat4<T>::det() const {

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
m::tmat4<T> m::tmat4<T>::cofactors() const {

    tmat4<T> result;
    T s = 1;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = s * getMinor(x, y).det();
            s *= -1;
        }

        s *= -1 * ((N + 1) % 2);
    }

    return result;
}

template <typename T>
m::tmat4<T> m::tmat4<T>::transpose() const {

    tmat4<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = get(y, x);
        }
    }

    return result;
}

template <typename T>
m::tmat4<T> m::tmat4<T>::adjoint() const {

    return cofactors().transpose();
}

template <typename T>
m::tmat4<T> m::tmat4<T>::inverse() const {

    T determinant = det();

    if (util::checkZero(determinant)) throw std::invalid_argument("Singular matrix can't be inverted");

    tmat4<T> adj = adjoint();

    return adj / determinant;
}

template <typename T>
m::tmat4<T> m::tmat4<T>::identity() {

    return tmat4<T>();
}

template <typename T>
m::tmat4<T> m::tmat4<T>::scale(T x, T y, T z) {

    return tmat4<T>{x, 0, 0, 0,
                    0, y, 0, 0,
                    0, 0, z, 0,
                    0, 0, 0, 1};
}

template <typename T>
m::tmat4<T> m::tmat4<T>::scale(m::tvec3<T> v) {

    return scale(v.x, v.y, v.z);
}

template <typename T>
m::tmat4<T> m::tmat4<T>::translate(T x, T y, T z) {

    return tmat4<T>{1, 0, 0, 0, // note: this layout is transposed because column-major
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    x, y, z, 1};
}

template <typename T>
m::tmat4<T> m::tmat4<T>::orthogonalProjection(T left, T right, T bottom, T top, T near, T far) { // note: eye space looks into -z, this projection flips that

    T rml = right - left;
    T rpl = right + left;
    T tmb = top - bottom;
    T tpb = top + bottom;
    T fmn = far - near;
    T fpn = far + near;

    return tmat4<T>::translate(-rpl/rml, -tpb/tmb, -fpn,fmn)
         * tmat4<T>::scale(2/rml, 2/tmb, -2/fmn);
}

template <typename T>
m::tmat4<T> m::tmat4<T>::perspectiveProjection(T fovHorizontal, T fovVertical, T near, T far) {

    T cotH = static_cast<T>(1.0 / std::tan(static_cast<double>(fovHorizontal) / 2.0));
    T cotV = static_cast<T>(1.0 / std::tan(static_cast<double>(fovVertical) / 2.0));

    return tmat4<T>{cotH, 0, 0, 0, // note: this layout is transposed because column-major
                    0, cotV, 0, 0,
                    0, 0, far/(far - near), far*near/(near - far),
                    0, 0, 1, 0};
}

template <typename T>
m::tmat4<T> m::operator+(const m::tmat4<T> &a, const m::tmat4<T> &b) {

    tmat4<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] += b.values[i];

    return result;
}

template <typename T>
m::tmat4<T> m::operator-(const tmat4<T> &a, const m::tmat4<T> &b) {

    tmat4<T> result(a);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] -= b.values[i];

    return result;
}

template <typename T>
m::tmat4<T> m::operator*(T scalar, const m::tmat4<T> &matrix) {

    tmat4<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] *= scalar;

    return result;
}

template <typename T>
m::tmat4<T> m::operator*(const m::tmat4<T> &matrix, T scalar) {

    return scalar * matrix;
}

template <typename T>
m::tmat4<T> m::operator/(const m::tmat4<T> &matrix, T scalar) {

    tmat4<T> result(matrix);

    for (size_t i = 0; i < N * N; i++)
        result.values[i] /= scalar;

    return result;
}

template <typename T>
m::tmat4<T> m::operator*(const m::tmat4<T> &a, const tmat4<T> &b) {

    tmat4<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = tvec4<T>::dot(a.getRow(y), b.getColumn(x));
        }
    }

    return result;
}

template <typename T>
m::tmat4<T> m::operator/(const m::tmat4<T> &a, const m::tmat4<T> &b) {

    return b.inverse() * a;
}

template <typename T>
m::tvec4<T> m::operator*(const m::tmat4<T> &matrix, const m::tvec4<T> &vector) {

    tvec4<T> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(y) = tvec4<T>::dot(matrix.getRow(y), vector);
        }
    }

    return result;
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const m::tmat4<T> &matrix) {

    return stream << std::endl << "__\t\t\t\t\t\b \b__" << std::endl
                               << "|\t" << matrix.values[0] << "\t" << matrix.values[4] << "\t" << matrix.values[8]  << "\t" << matrix.values[12] << "\t|" << std::endl
                               << "|\t" << matrix.values[1] << "\t" << matrix.values[5] << "\t" << matrix.values[9]  << "\t" << matrix.values[13] << "\t|" << std::endl
                               << "|\t" << matrix.values[2] << "\t" << matrix.values[6] << "\t" << matrix.values[10] << "\t" << matrix.values[14] << "\t|" << std::endl
                               << "|\t" << matrix.values[3] << "\t" << matrix.values[7] << "\t" << matrix.values[11] << "\t" << matrix.values[15] << "\t|" << std::endl
                               << "--\t\t\t\t\t\b \b--";
}

#undef N



#endif
