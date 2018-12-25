#ifndef __m_quat_h__
#define __m_quat_h__

#include <cmath>
#include <ostream>

#ifndef __m_vec_h__
#include <m/vec.h>
#endif

namespace m {

    template <typename T>
    class tquat {

    private:

        T w;
        tvec<T, 3> ijk;

    public:

#define tvec3 tvec<T, 3>

#define BINDING(name, value, type) const type & name () const { return value ; } \
                                         type & name ()       { return value ; }

        BINDING(real, w, T)
        BINDING(imaginary, ijk, tvec3)
        BINDING(i, ijk.x(), T)
        BINDING(j, ijk.y(), T)
        BINDING(k, ijk.z(), T)

#undef BINDING

#undef tvec3

        tquat() : w(0), ijk(0, 0, 0) {}

        tquat(T a, T b, T c, T d) : w(a), ijk(b, c, d) {}

        tquat(const tquat &other) : w(other.w), ijk(other.ijk) {}

        tquat(T r) : w(r), ijk(0, 0, 0) {}

        tquat(const tvec<T, 3> &vector)  :w(0), ijk(vector) {}

        tquat(T r, const tvec<T, 3> &vector) : w(r), ijk(vector) {}

        T magnSqr() const {

            return w * w + ijk.magnSqr();
        }

        double magn() const {

            return std::sqrt(static_cast<double>(magnSqr()));
        }

        tquat<T> conjugate() const {

            return tquat<T>(w, -ijk);
        }

        tquat<T> inverse() const {

            T ls = magnSqr();

            if (util::checkZero(ls)) throw std::invalid_argument("zero quaternion has no inverse");

            return conjugate() / ls;
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
                                                                 //                        /-
        static tquat rotation(T angle, const tvec<T, 3> &axis) { // NOTE: Right-handed: ---|--> axis
                                                                 //                        \-> rotation
            double halfAngle = static_cast<double>(angle) / 2.0;

            T c = static_cast<T>(std::cos(halfAngle));
            T s = static_cast<T>(std::sin(halfAngle));

            return c + s * tquat<T>(axis);
        }

        friend tquat<T> operator+(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.w + b.w, a.ijk + b.ijk);
        }

        friend tquat<T> operator+(const tquat<T> &quaternion, T scalar) {

            return tquat<T>(quaternion.w + scalar, quaternion.ijk);
        }

        friend tquat<T> operator+(T scalar, const tquat<T> &quaternion) {

            return quaternion + scalar;
        }

        friend tquat<T> operator-(const tquat<T> &a) {

            return tquat<T>(-a.w, -a.ijk);
        }

        friend tquat<T> operator-(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.w - b.w, a.ijk - b.ijk);
        }

        friend tquat<T> operator-(const tquat<T> &quaternion, T scalar) {

            return tquat<T>(quaternion.w - scalar, quaternion.ijk);
        }

        friend tquat<T> operator-(T scalar, const tquat<T> &quaternion) {

            return tquat<T>(scalar - quaternion.w, -quaternion.ijk);
        }

        friend tquat<T> operator*(const tquat<T> &a, const tquat<T> &b) {

            return tquat<T>(a.real() * b.i() + a.i() * b.real() + a.j() * b.k() - a.k() * b.j(),
                            a.real() * b.j() - a.i() * b.k() + a.j() * b.real() + a.k() * b.i(),
                            a.real() * b.k() + a.i() * b.j() - a.j() * b.i() + a.k() * b.real(),
                            a.real() * b.real() - a.i() * b.i() - a.j() * b.j() - a.k() * b.k());
        }

        friend tquat<T> operator*(T scalar, const tquat<T> &quaternion) {

            return tquat<T>(scalar * quaternion.w, scalar * quaternion.ijk);
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

            return tquat<T>(quaternion.w / scalar, quaternion.ijk / scalar);
        }

        friend std::ostream &operator<<(std::ostream &stream, const tquat<T> &quaternion) {

            bool nonZero = false;

            stream << std::fixed << std::setprecision(

#ifdef m_PRECISION

            m_PRECISION

#else
                    
            2

#endif
                    
            ) << "(";

            if (!util::checkZero(quaternion.real())) {

                stream << quaternion.real();
                nonZero = true;
            }

            if (!util::checkZero(quaternion.i())) {

                if (nonZero) {

                    stream << (quaternion.i() > 0 ? " + " : " - ");
                    stream << std::abs(quaternion.i());

                } else {

                    stream << quaternion.i();
                    nonZero = true;
                }

                stream << "i";
            }

            if (!util::checkZero(quaternion.j())) {
                
                if (nonZero) {

                    stream << (quaternion.j() > 0 ? " + " : " - ");
                    stream << std::abs(quaternion.j());

                } else {

                    stream << quaternion.j();
                    nonZero = true;
                }

                stream << "j";
            }

            if (!util::checkZero(quaternion.k())) {
                
                if (nonZero) {

                    stream << (quaternion.k() > 0 ? " + " : " - ");
                    stream << std::abs(quaternion.k());

                } else {

                    stream << quaternion.k();
                }

                stream << "k";
            }

            if (!nonZero) stream << "0";

            return stream << ")";
        }

        tquat<T> &operator+=(const tquat<T> &other) {

            *this = *this + other;

            return *this;
        }

        tquat<T> &operator+=(T other) {

            *this = *this + other;

            return *this;
        }

        tquat<T> &operator-=(const tquat<T> &other) {

            *this = *this - other;

            return *this;
        }

        tquat<T> &operator-=(T other) {

            *this = *this - other;

            return *this;
        }

        tquat<T> &operator*=(const tquat<T> &other) {

            *this = *this * other;

            return *this;
        }

        tquat<T> &operator*=(T other) {

            *this = *this * other;

            return *this;
        }

        tquat<T> &operator/=(const tquat<T> &other) {

            *this = *this / other;

            return *this;
        }

        tquat<T> &operator/=(T other) {

            *this = *this / other;

            return *this;
        }
    };

    typedef tquat<int>    iquat;
    typedef tquat<long>   lquat;
    typedef tquat<float>  quat;
    typedef tquat<double> dquat;
}

#endif
