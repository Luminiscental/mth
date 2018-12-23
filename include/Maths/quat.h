#ifndef __Maths_tquat_h__
#define __Maths_tquat_h__

#include <cmath>
#include <ostream>

#ifndef __Maths_vec_h__
#include <Maths/vec.h>
#endif

namespace m {

    template <typename T>
    struct tquat {

        T x, y, z, w; // w + ix + jy + kz

        T real() const {
            
            return w; 
        }

        T i() const {
            
            return x; 
        }

        T j() const {
            
            return y; 
        }

        T k() const {
            
            return z; 
        }

        tvec<T, 3> imaginary() const {
            
            return tvec<T, 3>(i(), j(), k()); 
        }

        tquat() : x(0), y(0), z(0), w(1) {}
        tquat(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {}
        tquat(const tquat &other) : x(other.x), y(other.y), z(other.z), w(other.w) {}
        tquat(T real) : x(0), y(0), z(0), w(real) {}
        tquat(const tvec<T, 3> &imaginary) : x(imaginary.x), y(imaginary.y), z(imaginary.z), w(0) {}

        T magnSqr() const {

            return x * x + y * y + z * z + w * w;
        }

        double magn() const {

            return std::sqrt(static_cast<double>(magnSqr()));
        }

        tquat<T> conjugate() const {

            return tquat<T>(-x, -y, -z, w);
        }

        tquat<T> inverse() const {

            return conjugate() / magnSqr();
        }

        tquat<T> unit() const {

            T l = static_cast<T>(magn());

            if (util::checkZero(l)) throw std::invalid_argument("zero quaternion has no unit equivalent");

            return *this / l;
        }

        tvec<T, 3> rotate(const m::tvec<T, 3> &vector) const {

            return (*this * tquat<T>(vector) * inverse()).imaginary();
        }

        static tquat identity() {
            
            return tquat(); 
        }

        static tquat rotation(T angle, const tvec<T, 3> &axis) {

            double halfAngle = static_cast<double>(angle) / 2.0;
            return static_cast<T>(std::cos(halfAngle)) + static_cast<T>(std::sin(halfAngle)) * tquat<T>(axis);
        }

        friend tquat<T> operator+(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
        }

        friend tquat<T> operator+(const tquat<T> &quaternion, T scalar) {

            return tquat<T>(scalar) + quaternion;
        }

        friend tquat<T> operator+(T scalar, const tquat<T> &quaternion) {

            return quaternion + scalar;
        }

        friend tquat<T> operator-(const tquat<T> &a) {

            return tquat<T>(-a.x, -a.y, -a.z, -a.w);
        }

        friend tquat<T> operator-(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
        }

        friend tquat<T> operator-(const tquat<T> &quaternion, T scalar) {

            return quaternion - tquat<T>(scalar);
        }

        friend tquat<T> operator-(T scalar, const tquat<T> &quaternion) {

            return tquat<T>(scalar) - quaternion;
        }

        friend tquat<T> operator*(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
                            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
                            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
                            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
        }

        friend tquat<T> operator*(T scalar, const tquat<T> &quaternion) {

            return tquat<T>(scalar) * quaternion;
        }

        friend tquat<T> operator*(const tquat<T> &quaternion, T scalar) {

            return scalar * quaternion;
        }

        friend tquat<T> operator/(const tquat<T> &a, const tquat<T> &b) {

            return a * b.inverse();
        }

        friend tquat<T> operator/(T scalar, const tquat<T> &quaternion) {

            return tquat<T>(scalar) / quaternion;
        }

        friend tquat<T> operator/(const tquat<T> &quaternion, T scalar) {

            return quaternion / tquat<T>(scalar);
        }

        friend std::ostream &operator<<(std::ostream &stream, const tquat<T> &quaternion) {

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
    };

    typedef tquat<int>    iquat;
    typedef tquat<long>   lquat;
    typedef tquat<float>  quat;
    typedef tquat<double> dquat;
}

#endif
