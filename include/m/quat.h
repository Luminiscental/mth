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

#define BINDING(name, value) const auto & name () const { return value ; } \
                             auto & name () { return value ; }

        BINDING(real, w)
        BINDING(imaginary, ijk)
        BINDING(i, ijk.x())
        BINDING(j, ijk.y())
        BINDING(k, ijk.z())

#undef BINDING

        constexpr tquat() noexcept
            :w(0), ijk(0, 0, 0) {}

        constexpr tquat(T a, T b, T c, T d) noexcept
            :w(a), ijk(b, c, d) {}

        constexpr tquat(T r) noexcept
            :w(r), ijk(0, 0, 0) {}

        constexpr tquat(const tvec<T, 3> &vector) noexcept
            :w(0), ijk(vector) {}

        constexpr tquat(T r, const tvec<T, 3> &vector) noexcept
            :w(r), ijk(vector) {}

        auto magnSqr() const noexcept {

            return w * w + ijk.magnSqr();
        }

        auto magn() const noexcept {

            return std::sqrt(static_cast<double>(magnSqr()));
        }

        constexpr auto conjugate() const noexcept {

            return tquat<T>(w, -ijk);
        }

        auto inverse() const {

            auto ls = magnSqr();

            if (util::checkZero(ls)) throw std::invalid_argument("zero quaternion has no inverse");

            return conjugate() / ls;
        }

        auto unit() const {

            auto l = static_cast<T>(magn());

            if (util::checkZero(l)) throw std::invalid_argument("zero quaternion has no unit equivalent");

            return *this / l;
        }

        auto rotate(const m::tvec<T, 3> &vector) const {

            return (*this * tquat<T>(vector) * inverse()).imaginary();
        }

        static auto identity() {
            
            return tquat(); 
        }
                                                                //                        /-
        static auto rotation(T angle, const tvec<T, 3> &axis) { // NOTE: Right-handed: ---|--> axis
                                                                //                        \-> rotation
            auto halfAngle = static_cast<double>(angle) / 2.0;

            auto c = static_cast<T>(std::cos(halfAngle));
            auto s = static_cast<T>(std::sin(halfAngle));

            return c + s * tquat<T>(axis);
        }

        auto &operator+=(const tquat<T> &rhs) {

            w += rhs.w;
            ijk += rhs.ijk;

            return *this;
        }

        auto &operator+=(T rhs) {

            w += rhs;

            return *this;
        }

        auto &operator-=(const tquat<T> &rhs) {

            w -= rhs.w;
            ijk -= rhs.ijk;

            return *this;
        }

        auto &operator-=(T rhs) {

            w -= rhs;

            return *this;
        }

        auto &operator*=(const tquat<T> &rhs) {

            *this = *this * rhs;

            return *this;
        }

        auto &operator*=(T rhs) {

            w *= rhs;
            ijk *= rhs;

            return *this;
        }

        auto &operator/=(const tquat<T> &rhs) {

            *this = *this / rhs;

            return *this;
        }

        auto &operator/=(T rhs) {

            w /= rhs;
            ijk /= rhs;

            return *this;
        }

        friend auto operator+(const tquat<T> &lhs, const tquat<T> &rhs) {

            tquat<T> result = lhs;

            return result += rhs;
        }

        friend auto operator+(const tquat<T> &lhs, T rhs) {

            tquat<T> result = lhs;

            return result += rhs;
        }

        friend auto operator+(T lhs, const tquat<T> &rhs) {

            return rhs + lhs;
        }

        friend auto operator-(const tquat<T> &rhs) {

            return tquat<T>(-rhs.w, -rhs.ijk);
        }

        friend auto operator-(const tquat<T> &lhs, const tquat<T> &rhs) {

            tquat<T> result = lhs;

            return result -= rhs;
        }

        friend auto operator-(const tquat<T> &lhs, T rhs) {

            tquat<T> result = lhs;

            return result -= rhs;
        }

        friend auto operator-(T lhs, const tquat<T> &rhs) {

            tquat<T> result = lhs;

            return result -= rhs;
        }

        friend auto operator*(const tquat<T> &lhs, const tquat<T> &rhs) {

            return tquat<T>(lhs.real() * rhs.i() + lhs.i() * rhs.real() + lhs.j() * rhs.k() - lhs.k() * rhs.j(),
                            lhs.real() * rhs.j() - lhs.i() * rhs.k() + lhs.j() * rhs.real() + lhs.k() * rhs.i(),
                            lhs.real() * rhs.k() + lhs.i() * rhs.j() - lhs.j() * rhs.i() + lhs.k() * rhs.real(),
                            lhs.real() * rhs.real() - lhs.i() * rhs.i() - lhs.j() * rhs.j() - lhs.k() * rhs.k());
        }

        friend auto operator*(T lhs, const tquat<T> &rhs) {

            tquat<T> result = rhs;

            return result *= lhs;
        }

        friend auto operator*(const tquat<T> &lhs, T rhs) {

            return rhs * lhs;
        }

        friend auto operator/(const tquat<T> &lhs, const tquat<T> &rhs) {

            return lhs * rhs.inverse();
        }

        friend auto operator/(T lhs, const tquat<T> &rhs) {

            return tquat<T>(lhs) / rhs;
        }

        friend auto operator/(const tquat<T> &lhs, T rhs) {

            tquat<T> result = lhs;

            return result /= rhs;
        }

        friend auto &operator<<(std::ostream &lhs, const tquat<T> &rhs) {

            bool nonZero = false;

            lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

            if (!util::checkZero(rhs.real())) {

                lhs << rhs.real();
                nonZero = true;
            }

            if (!util::checkZero(rhs.i())) {

                if (nonZero) {

                    lhs << (rhs.i() > 0 ? " + " : " - ");
                    lhs << std::abs(rhs.i());

                } else {

                    lhs << rhs.i();
                    nonZero = true;
                }

                lhs << "i";
            }

            if (!util::checkZero(rhs.j())) {
                
                if (nonZero) {

                    lhs << (rhs.j() > 0 ? " + " : " - ");
                    lhs << std::abs(rhs.j());

                } else {

                    lhs << rhs.j();
                    nonZero = true;
                }

                lhs << "j";
            }

            if (!util::checkZero(rhs.k())) {
                
                if (nonZero) {

                    lhs << (rhs.k() > 0 ? " + " : " - ");
                    lhs << std::abs(rhs.k());

                } else {

                    lhs << rhs.k();
                }

                lhs << "k";
            }

            if (!nonZero) lhs << "0";

            return lhs << ")";
        }
    };

    using iquat = tquat<int>;
    using lquat = tquat<long>;
    using  quat = tquat<float>;
    using dquat = tquat<double>;
}

#endif
