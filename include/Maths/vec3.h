#ifndef __Maths_vec3_h__
#define __Maths_vec3_h__

#include <cmath>
#include <ostream>

#ifndef __Maths_vec2_h__
#include <Maths/vec2.h>
#endif

namespace m {

    template <typename T> // requires T + T, T - T, -T, T * T, T / T, static_cast<double>(T), std::ostream << T
    struct tvec3 {

        T x, y, z;

        tvec3() :x(0), y(0), z(0) {}
        tvec3(T x, T y, T z) :x(x), y(y), z(z) {}
        tvec3(const tvec3 &other) :x(other.x), y(other.y), z(other.z) {}
        tvec3(const tvec2<T> &xy, T z) :x(xy.x), y(xy.y), z(z) {}
        tvec3(T x, const tvec2<T> &yz) :x(x), y(yz.x), z(yz.y) {}

        tvec2<T> xy() const;
        tvec2<T> yz() const;

        T magnSqr() const;
        double magn() const;

        static T dot(const tvec3 &a, const tvec3 &b);
        static tvec3 cross(const tvec3 &a, const tvec3 &b);
    };

    typedef tvec3<int>    ivec3;
    typedef tvec3<long>   lvec3;
    typedef tvec3<float>   vec3;
    typedef tvec3<double> dvec3;

    template <typename T>
    bool operator==(const tvec3<T> &a, const tvec3<T> &b);

    template <typename T>
    bool operator!=(const tvec3<T> &a, const tvec3<T> &b);

    template <typename T>
    tvec3<T> operator+(const tvec3<T> &a, const tvec3<T> &b);

    template <typename T>
    tvec3<T> operator-(const tvec3<T> &a, const tvec3<T> &b);

    template <typename T>
    tvec3<T> operator-(const tvec3<T> &vector);

    template <typename T>
    tvec3<T> operator*(T scalar, const tvec3<T> &vector);

    template <typename T>
    tvec3<T> operator*(const tvec3<T> &vector, T scalar);

    template <typename T>
    tvec3<T> operator/(const tvec3<T> &vector, T scalar);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tvec3<T> &vector);
}

// Template implementation

template <typename T>
m::tvec2<T> m::tvec3<T>::xy() const {

    return m::tvec2<T>(x, y);
}

template <typename T>
m::tvec2<T> m::tvec3<T>::yz() const {

    return m::tvec2<T>(y, z);
}

template <typename T>
T m::tvec3<T>::magnSqr() const {

    return x * x + y * y + z * z;
}

template <typename T>
double m::tvec3<T>::magn() const {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
T m::tvec3<T>::dot(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template <typename T>
m::tvec3<T> m::tvec3<T>::cross(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return m::tvec3<T>(a.y * b.z - a.z * b.y, a.x * b.z - a.z * b.x, a.x * b.y - a.y * b.x);
}

template <typename T>
bool m::operator==(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return a.xy() == b.xy() && a.z == b.z;
}

template <typename T>
bool m::operator!=(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return a.xy() != b.xy() || a.z != b.z;
}

template <typename T>
m::tvec3<T> m::operator+(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return m::tvec3<T>(a.xy() + b.xy(), a.z + b.z);
}

template <typename T>
m::tvec3<T> m::operator-(const m::tvec3<T> &a, const m::tvec3<T> &b) {

    return m::tvec3<T>(a.xy() - b.xy(), a.z - b.z);
}

template <typename T>
m::tvec3<T> m::operator-(const m::tvec3<T> &vector) {

    return m::tvec3<T>(-vector.xy(), -vector.z);
}

template <typename T>
m::tvec3<T> m::operator*(T scalar, const m::tvec3<T> &vector) {

    return m::tvec3<T>(scalar * vector.xy(), scalar * vector.z);
}

template <typename T>
m::tvec3<T> m::operator*(const m::tvec3<T> &vector, T scalar) {

    return scalar * vector;
}

template <typename T>
m::tvec3<T> m::operator/(const m::tvec3<T> &vector, T scalar) {

    return m::tvec3<T>(vector.xy() / scalar, vector.z / scalar);
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const m::tvec3<T> &vector) {

    return stream << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
}

#endif
