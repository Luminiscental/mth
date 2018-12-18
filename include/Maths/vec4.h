#ifndef __Maths_vec4_h__
#define __Maths_vec4_h__

#include <cmath>
#include <ostream>
#include <functional>

#ifndef __Maths_vec3_h__
#include <Maths/vec3.h>
#endif

#ifndef __Maths_vec2_h__
#include <Maths/vec2.h>
#endif

namespace m {

    template <typename T> // requires 0, T + T, T - T, -T, T * T, T / T, static_cast<double>(T), std::ostream << T
    struct tvec4 {

        T x, y, z, w;

        tvec4() :x(0), y(0), z(0), w(0) {}
        tvec4(T x, T y, T z, T w) :x(x), y(y), z(z), w(w) {}
        tvec4(const tvec4 &other) :x(other.x), y(other.y), z(other.z), w(other.w) {}
        tvec4(const tvec3<T> &xyz, T w) :x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
        tvec4(T x, const tvec3<T> &yzw) :x(x), y(yzw.x), z(yzw.y), w(yzw.z) {}
        tvec4(const tvec2<T> &xy, const tvec2<T> &zw) :x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
        tvec4(const tvec2<T> &xy, T z, T w) :x(xy.x), y(xy.y), z(z), w(w) {}
        tvec4(T x, T y, const tvec2<T> &zw) :x(x), y(y), z(zw.x), w(zw.y) {}
        tvec4(T x, const tvec2<T> &yz, T w) :x(x), y(yz.x), z(yz.y), w(w) {}

        tvec2<T> xy() const;
        tvec2<T> yz() const;
        tvec2<T> zw() const;

        tvec3<T> xyz() const;
        tvec3<T> yzw() const;

        void forEach(std::function<void(const T&)> action) const;
        void forEach(std::function<void(T&)> action);

        T magnSqr() const;
        double magn() const;

        static T dot(const tvec4 &a, const tvec4 &b);
    };

    typedef tvec4<int>    ivec4;
    typedef tvec4<long>   lvec4;
    typedef tvec4<float>  vec4;
    typedef tvec4<double> dvec4;

    template <typename T>
    bool operator==(const tvec4<T> &a, const tvec4<T> &b);

    template <typename T>
    bool operator!=(const tvec4<T> &a, const tvec4<T> &b);

    template <typename T>
    tvec4<T> operator+(const tvec4<T> &a, const tvec4<T> &b);

    template <typename T>
    tvec4<T> operator-(const tvec4<T> &a, const tvec4<T> &b);

    template <typename T>
    tvec4<T> operator-(const tvec4<T> &vector);

    template <typename T>
    tvec4<T> operator*(T scalar, const tvec4<T> &vector);

    template <typename T>
    tvec4<T> operator*(const tvec4<T> &vector, T scalar);

    template <typename T>
    tvec4<T> operator/(const tvec4<T> &vector, T scalar);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tvec4<T> &vector);
}

// Template implementation

template <typename T>
void m::tvec4<T>::forEach(std::function<void(const T&)> action) const { 

    xyz().forEach(action);
    action(w);
}

template <typename T>
void m::tvec4<T>::forEach(std::function<void(T&)> action) {

    xyz().forEach(action);
    action(w);
}

template <typename T>
m::tvec2<T> m::tvec4<T>::xy() const {

    return m::tvec2<T>(x, y);
}

template <typename T>
m::tvec2<T> m::tvec4<T>::yz() const {

    return m::tvec2<T>(y, z);
}

template <typename T>
m::tvec2<T> m::tvec4<T>::zw() const {

    return m::tvec2<T>(z, w);
}

template <typename T>
m::tvec3<T> m::tvec4<T>::xyz() const {

    return m::tvec3<T>(x, y, z);
}

template <typename T>
m::tvec3<T> m::tvec4<T>::yzw() const {

    return m::tvec3<T>(y, z, w);
}

template <typename T>
T m::tvec4<T>::magnSqr() const {

    return x * x + y * y + z * z + w * w;
}

template <typename T>
double m::tvec4<T>::magn() const {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
T m::tvec4<T>::dot(const m::tvec4<T> &a, const m::tvec4<T> &b) {

    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

template <typename T>
bool m::operator==(const m::tvec4<T> &a, const m::tvec4<T> &b) {

    return a.xyz() == b.xyz() && a.w == b.w;
}

template <typename T>
bool m::operator!=(const m::tvec4<T> &a, const m::tvec4<T> &b) {

    return a.xyz() != b.xyz() || a.w != b.w;
}

template <typename T>
m::tvec4<T> m::operator+(const m::tvec4<T> &a, const m::tvec4<T> &b) {

    return m::tvec4<T>(a.xyz() + b.xyz(), a.w + b.w);
}

template <typename T>
m::tvec4<T> m::operator-(const m::tvec4<T> &a, const m::tvec4<T> &b) {

    return m::tvec4<T>(a.xyz() - b.xyz(), a.w - b.w);
}

template <typename T>
m::tvec4<T> m::operator-(const m::tvec4<T> &vector) {

    return m::tvec4<T>(-vector.xyz(), -vector.w);
}

template <typename T>
m::tvec4<T> m::operator*(T scalar, const m::tvec4<T> &vector) {

    return m::tvec4<T>(scalar * vector.xyz(), scalar * vector.w);
}

template <typename T>
m::tvec4<T> m::operator*(const m::tvec4<T> &vector, T scalar) {

    return scalar * vector;
}

template <typename T>
m::tvec4<T> m::operator/(const m::tvec4<T> &vector, T scalar) {

    return m::tvec4<T>(vector.xyz() / scalar, vector.w / scalar);
}

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const tvec4<T> &vector) {

    return stream << "(" << vector.x << ", " << vector.y << ", " << vector.z << ", " << vector.w << ")";
}

#endif
