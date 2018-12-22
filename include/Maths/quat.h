#ifndef __Maths_tquat_h__
#define __Maths_tquat_h__

#include <cmath>
#include <ostream>

#ifndef __Maths_tvec3_h__
#include <Maths/vec3.h>
#endif

namespace m {

    // requires 0, 1, T + T, T - T, -T, T * T, T / T, static_cast<double>(T), static_cast<T>(double), std::ostream << T, std::abs(T), std::numeric_limits<T>::epsilon()
    template <typename T>
    struct tquat {

        T x, y, z, w; // w + i.x + j.y + k.z

        inline T real() const { return w; }
        inline T i() const { return x; }
        inline T j() const { return y; }
        inline T k() const { return z; }
        inline tvec3<T> imaginary() const { return tvec3<T>(i(), j(), k()); }

        tquat() :x(0), y(0), z(0), w(1) {}
        tquat(T x, T y, T z, T w) :x(x), y(y), z(z), w(w) {}
        tquat(const tquat &other) :x(other.x), y(other.y), z(other.z), w(other.w) {}
        tquat(T real) :x(0), y(0), z(0), w(real) {}
        tquat(const tvec3<T> &imaginary) :x(imaginary.x), y(imaginary.y), z(imaginary.z), w(0) {}

        T magnSqr() const;
        double magn() const;

        tquat conjugate() const;
        tquat inverse() const;
        tquat unit() const;

        tvec3<T> rotate(const tvec3<T> &vector) const;

        inline static tquat identity() { return tquat(); }
        static tquat rotation(T angle, const tvec3<T> &axis);
    };

    typedef tquat<int>    iquat;
    typedef tquat<long>   lquat;
    typedef tquat<float>  quat;
    typedef tquat<double> dquat;

    template <typename T>
    tquat<T> operator+(const tquat<T> &a, const tquat<T> &b);

    template <typename T>
    tquat<T> operator+(T scalar, const tquat<T> &quaternion);

    template <typename T>
    tquat<T> operator+(const tquat<T> &quaternion, T scalar);

    template <typename T>
    tquat<T> operator-(const tquat<T> &a);

    template <typename T>
    tquat<T> operator-(const tquat<T> &a, const tquat<T> &b);

    template <typename T>
    tquat<T> operator-(const tquat<T> &quaternion, T scalar);

    template <typename T>
    tquat<T> operator-(T scalar, const tquat<T> &quaternion);

    template <typename T>
    tquat<T> operator*(const tquat<T> &a, const tquat<T> &b);

    template <typename T>
    tquat<T> operator*(T scalar, const tquat<T> &quaternion);

    template <typename T>
    tquat<T> operator*(const tquat<T> &quaternion, T scalar);

    template <typename T>
    tquat<T> operator/(const tquat<T> &a, const tquat<T> &b);

    template <typename T>
    tquat<T> operator/(const tquat<T> &quaternion, T scalar);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tquat<T> &quaternion);
}

// Template implementation

template <typename T>
T m::tquat<T>::magnSqr() const {

    return x * x + y * y + z * z + w * w;
}

template <typename T>
double m::tquat<T>::magn() const {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
m::tquat<T> m::tquat<T>::conjugate() const {

    return tquat<T>(-x, -y, -z, w);
}

template <typename T>
m::tquat<T> m::tquat<T>::inverse() const {

    return conjugate() / magnSqr();
}

template <typename T>
m::tquat<T> m::tquat<T>::unit() const {

    return *this / static_cast<T>(magn());
}

template <typename T>
m::tvec3<T> m::tquat<T>::rotate(const m::tvec3<T> &vector) const {

    return (*this * tquat<T>(vector) * inverse()).imaginary();
}

template <typename T>
m::tquat<T> m::tquat<T>::rotation(T angle, const m::tvec3<T> &axis) {

    double halfAngle = static_cast<double>(angle) / 2.0;
    return static_cast<T>(std::cos(halfAngle)) + static_cast<T>(std::sin(halfAngle)) * tquat<T>(axis);
}

template <typename T>
m::tquat<T> m::operator+(const m::tquat<T> &a, const m::tquat<T> &b) {

    return tquat<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

template <typename T>
m::tquat<T> m::operator+(const m::tquat<T> &quaternion, T scalar) {

    return tquat<T>(scalar) + quaternion;
}

template <typename T>
m::tquat<T> m::operator+(T scalar, const m::tquat<T> &quaternion) {

    return quaternion + scalar;
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &a) {

    return tquat<T>(-a.x, -a.y, -a.z, -a.w);
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &a, const m::tquat<T> &b) {

    return tquat<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &quaternion, T scalar) {

    return quaternion - tquat<T>(scalar);
}

template <typename T>
m::tquat<T> m::operator-(T scalar, const m::tquat<T> &quaternion) {

    return tquat<T>(scalar) - quaternion;
}

template <typename T>
m::tquat<T> m::operator*(const m::tquat<T> &a, const m::tquat<T> &b) {

    return tquat<T>(a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
                    a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
                    a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
                    a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
}

template <typename T>
m::tquat<T> m::operator*(T scalar, const m::tquat<T> &quaternion) {

    return tquat<T>(scalar * quaternion.x, scalar * quaternion.y, scalar * quaternion.z, scalar * quaternion.w);
}

template <typename T>
m::tquat<T> m::operator*(const m::tquat<T> &quaternion, T scalar) {

    return scalar * quaternion;
}

template <typename T>
m::tquat<T> m::operator/(const m::tquat<T> &a, const m::tquat<T> &b) {

    return a * b.inverse();
}

template <typename T>
m::tquat<T> m::operator/(const m::tquat<T> &quaternion, T scalar) {

    return tquat<T>(quaternion.x / scalar, quaternion.y / scalar, quaternion.z / scalar, quaternion.w / scalar);
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const m::tquat<T> &quaternion) {

    bool nonZero = false;

    stream << "(";

    if (!util::checkZero(quaternion.w)) {

        stream << quaternion.w;
        nonZero = true;
    }

    if (!util::checkZero(quaternion.x)) {

        if (nonZero) stream << " + ";
        stream << quaternion.x << "i";
        nonZero = true;
    }

    if (!util::checkZero(quaternion.y)) {

        if (nonZero) stream << " + ";
        stream << quaternion.y << "j";
        nonZero = true;
    }

    if (!util::checkZero(quaternion.z)) {

        if (nonZero) stream << " + ";
        stream << quaternion.z << "k";
    }

    if (!nonZero) stream << "0";

    return stream << ")";
}

#endif
