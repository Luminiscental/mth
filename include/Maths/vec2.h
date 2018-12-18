#ifndef __Maths_vec2_h__
#define __Maths_vec2_h__

#include <cmath>
#include <ostream>

namespace m {

    template <typename T> // requires T + T, T - T, -T, T * T, T / T, static_cast<double>(T), std::ostream << T
    struct tvec2 {

        T x, y;

        tvec2() :x(0), y(0) {}
        tvec2(T x, T y) :x(x), y(y) {}
        tvec2(const tvec2 &other) :x(other.x), y(other.y) {}

        T magnSqr() const;
        double magn() const;

        double arg() const; // (-pi, pi]

        static T dot(const tvec2 &a, const tvec2 &b);
    };

    typedef tvec2<int>    ivec2;
    typedef tvec2<long>   lvec2;
    typedef tvec2<float>   vec2;
    typedef tvec2<double> dvec2;

    template <typename T>
    bool operator==(const tvec2<T> &a, const tvec2<T> &b);

    template <typename T>
    bool operator!=(const tvec2<T> &a, const tvec2<T> &b);

    template <typename T>
    tvec2<T> operator+(const tvec2<T> &a, const tvec2<T> &b);

    template <typename T>
    tvec2<T> operator-(const tvec2<T> &a, const tvec2<T> &b);

    template <typename T>
    tvec2<T> operator-(const tvec2<T> &a);

    template <typename T>
    tvec2<T> operator*(float scalar, const tvec2<T> &vector);

    template <typename T>
    tvec2<T> operator*(const tvec2<T> &vector, float scalar);

    template <typename T>
    tvec2<T> operator/(const tvec2<T> &vector, float scalar);

    template <typename T>
    std::ostream &operator<<(std::ostream &stream, const tvec2<T> &vector);
}

// Template implementation

template <typename T>
T m::tvec2<T>::magnSqr() const {

    return x * x + y * y;
}

template <typename T>
double m::tvec2<T>::magn() const {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
double m::tvec2<T>::arg() const {

    return std::atan2(y, x);
}

template <typename T>
T m::tvec2<T>::dot(const m::tvec2<T> &a, const m::tvec2<T> &b) {

    return a.x * b.x + a.y * b.y;
}

template <typename T>
bool m::operator==(const tvec2<T> &a, const tvec2<T> &b) {

    return a.x == b.x && a.y == b.y;
}

template <typename T>
bool m::operator!=(const tvec2<T> &a, const tvec2<T> &b) {

    return a.x != b.x || a.y != b.y;
}

template <typename T>
m::tvec2<T> m::operator+(const m::tvec2<T> &a, const m::tvec2<T> &b) {

    return m::tvec2<T>(a.x + b.x, a.y + b.y);
}

template <typename T>
m::tvec2<T> m::operator-(const m::tvec2<T> &a, const m::tvec2<T> &b) {

    return m::tvec2<T>(a.x - b.x, a.y - b.y);
}

template <typename T>
m::tvec2<T> m::operator-(const m::tvec2<T> &a) {

    return m::tvec2<T>(-a.x, -a.y);
}

template <typename T>
m::tvec2<T> m::operator*(float scalar, const m::tvec2<T> &vector) {

    return m::tvec2<T>(scalar * vector.x, scalar * vector.y);
}

template <typename T>
m::tvec2<T> m::operator*(const m::tvec2<T> &vector, float scalar) {

    return scalar * vector;
}

template <typename T>
m::tvec2<T> m::operator/(const m::tvec2<T> &vector, float scalar) {

    return m::tvec2<T>(vector.x / scalar, vector.y / scalar); }

template <typename T>
std::ostream &m::operator<<(std::ostream &stream, const tvec2<T> &vector) {

    return stream << "(" << vector.x << ", " << vector.y << ")";
}

#endif
