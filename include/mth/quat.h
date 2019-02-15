#ifndef mth_quat_h__
#define mth_quat_h__

/* <mth/quat.h> - quaternion header
 *      Includes the template class tquat representing a quaternion with coefficients of type T. Basic
 *      arithmetic operators are defined as well as member functions to get values such as modulus and
 *      reciprocal. Utility functions for creating and using rotation representations are also defined.
 */

#include <ostream>
#include <cmath>

#include <mth/mth.h>
#include <mth/vec.h>

namespace mth {

    template <typename T>
    class tquat {

    private:

        T w;
        tvec<T, 3> ijk;

    public:

#define BINDING(name, value)    constexpr const auto & name () const noexcept { return value; } \
                                constexpr       auto & name ()       noexcept { return value; }

        BINDING(real, w)
        BINDING(imaginary, ijk)
        BINDING(i, ijk.x())
        BINDING(j, ijk.y())
        BINDING(k, ijk.z())

#undef BINDING

        // Default initializes to zero
        constexpr tquat() noexcept 
            :w(0), ijk(0, 0, 0) {}

        // Order of coefficients is defined as real, i, j, k
        constexpr tquat(T a, T b, T c, T d) noexcept
            :w(a), ijk(b, c, d) {}

        // Initialize to a real number
        constexpr tquat(T r) noexcept
            :w(r), ijk(0, 0, 0) {}

        // Quaternion representation of a vector using ijk
        constexpr tquat(const tvec<T, 3> &vector) noexcept
            :w(0), ijk(vector) {}

        constexpr tquat(T r, const tvec<T, 3> &vector) noexcept
            :w(r), ijk(vector) {}

        template <typename U>
        constexpr operator tquat<U>() const noexcept {

            tquat<U> result;

            result.w = static_cast<U>(w);
            result.ijk = static_cast<tvec<U, 3>>(ijk);

            return result;
        }

        constexpr T magnSqr() const noexcept {

            return w * w + ijk.magnSqr();
        }

        // Converted to double for more accurate sqrt
        constexpr double magn() const noexcept {

            return std::sqrt(static_cast<double>(magnSqr()));
        }

        constexpr tquat<T> conjugate() const noexcept {

            return tquat<T>(w, -ijk);
        }

        // Returns the reciprocal
        constexpr tquat<T> inverse() const noexcept {

            auto ls = magnSqr();

            return conjugate() / ls;
        }

        constexpr tquat<T> unit() const noexcept {

            auto l = static_cast<T>(magn());

            return *this / l;
        }

        // Returns the rotation represented by *this applied to the given vector
        constexpr tvec<T, 3> rotate(const mth::tvec<T, 3> &vector) const noexcept {

            return (*this * tquat<T>(vector) * inverse()).imaginary();
        }

        // Returns 1
        constexpr static tquat<T> identity() noexcept {
            
            return tquat(); 
        }

        // Converts euler angle and axis to quaternion representation
        constexpr static tquat<T> rotation(T angle, const tvec<T, 3> &axis) noexcept {

            auto halfAngle = static_cast<double>(angle) / 2.0;

            auto c = static_cast<T>(std::cos(halfAngle));
            auto s = static_cast<T>(std::sin(halfAngle));

            return c + s * tquat<T>(axis);
        }

        constexpr tquat<T> &operator+=(const tquat<T> &rhs) noexcept {

            w += rhs.w;
            ijk += rhs.ijk;

            return *this;
        }

        constexpr tquat<T> &operator+=(T rhs) noexcept {

            w += rhs;

            return *this;
        }

        constexpr tquat<T> &operator-=(const tquat<T> &rhs) noexcept {

            w -= rhs.w;
            ijk -= rhs.ijk;

            return *this;
        }

        constexpr tquat<T> &operator-=(T rhs) noexcept {

            w -= rhs;

            return *this;
        }

        constexpr tquat<T> &operator*=(const tquat<T> &rhs) noexcept {

            *this = *this * rhs;

            return *this;
        }

        constexpr tquat<T> &operator*=(T rhs) noexcept {

            w *= rhs;
            ijk *= rhs;

            return *this;
        }

        constexpr tquat<T> &operator/=(const tquat<T> &rhs) noexcept {

            *this = *this / rhs;

            return *this;
        }

        constexpr tquat<T> &operator/=(T rhs) noexcept {

            w /= rhs;
            ijk /= rhs;

            return *this;
        }
    };

    template <typename T>
    constexpr tquat<T> operator+(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        tquat<T> result = lhs;

        return result += rhs;
    }

    template <typename T>
    constexpr tquat<T> operator+(const tquat<T> &lhs, T rhs) noexcept {

        tquat<T> result = lhs;

        return result += rhs;
    }

    template <typename T>
    constexpr tquat<T> operator+(T lhs, const tquat<T> &rhs) noexcept {

        return rhs + lhs;
    }

    template <typename T>
    constexpr tquat<T> operator-(const tquat<T> &rhs) noexcept {

        return tquat<T>(-rhs.w, -rhs.ijk);
    }

    template <typename T>
    constexpr tquat<T> operator-(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        tquat<T> result = lhs;

        return result -= rhs;
    }

    template <typename T>
    constexpr tquat<T> operator-(const tquat<T> &lhs, T rhs) noexcept {

        tquat<T> result = lhs;

        return result -= rhs;
    }

    template <typename T>
    constexpr tquat<T> operator-(T lhs, const tquat<T> &rhs) noexcept {

        tquat<T> result = lhs;

        return result -= rhs;
    }

    template <typename T>
    constexpr tquat<T> operator*(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        return tquat<T>(lhs.real() * rhs.real() - lhs.i() * rhs.i() - lhs.j() * rhs.j() - lhs.k() * rhs.k(),  // real
                        lhs.real() * rhs.i() + lhs.i() * rhs.real() + lhs.j() * rhs.k() - lhs.k() * rhs.j(),  // i
                        lhs.real() * rhs.j() - lhs.i() * rhs.k() + lhs.j() * rhs.real() + lhs.k() * rhs.i(),  // j
                        lhs.real() * rhs.k() + lhs.i() * rhs.j() - lhs.j() * rhs.i() + lhs.k() * rhs.real()); // k
    }

    template <typename T>
    constexpr tquat<T> operator*(T lhs, const tquat<T> &rhs) noexcept {

        tquat<T> result = rhs;

        return result *= lhs;
    }

    template <typename T>
    constexpr tquat<T> operator*(const tquat<T> &lhs, T rhs) noexcept {

        return rhs * lhs;
    }

    template <typename T>
    constexpr tquat<T> operator/(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        return lhs * rhs.inverse();
    }

    template <typename T>
    constexpr tquat<T> operator/(T lhs, const tquat<T> &rhs) noexcept {

        return tquat<T>(lhs) / rhs;
    }

    template <typename T>
    constexpr tquat<T> operator/(const tquat<T> &lhs, T rhs) noexcept {

        tquat<T> result = lhs;

        return result /= rhs;
    }

    template <typename T>
    constexpr bool operator==(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        return util::isEqual(lhs.real(), rhs.real()) && lhs.imag() == rhs.imag();
    }

    template <typename T>
    constexpr bool operator==(const T &lhs, const tquat<T> &rhs) noexcept {

        return util::isEqual(lhs, rhs.real()) && rhs.imag() == tvec<T, 3>(0, 0, 0);
    }

    template <typename T>
    constexpr bool operator==(const tquat<T> &lhs, const T &rhs) noexcept {

        return rhs == lhs;
    }

    template <typename T>
    constexpr bool operator!=(const tquat<T> &lhs, const tquat<T> &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T>
    constexpr bool operator!=(const T &lhs, const tquat<T> &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T>
    constexpr bool operator!=(const tquat<T> &lhs, const T &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &lhs, const tquat<T> &rhs) {

        bool realZero = util::isZero(rhs.real());
        bool iZero = util::isZero(rhs.i());
        bool jZero = util::isZero(rhs.j());
        bool kZero = util::isZero(rhs.k());

        bool zero = realZero && iZero && jZero && kZero;

        int termCount = !realZero + !iZero + !jZero + !kZero;
        bool multipleTerms = termCount > 1;

        bool termStreamed = false;

        if (zero) return lhs << "0";
        if (multipleTerms) lhs << "(";

        if (!realZero) {

            lhs << rhs.real();
            termStreamed = true;
        }

        if (!iZero) {

            if (termStreamed) {

                lhs << (rhs.i() > 0 ? " + " : " - ");
                lhs << std::abs(rhs.i());

            } else {

                lhs << rhs.i();
                termStreamed = true;
            }

            lhs << "i";
        }

        if (!jZero) {

            if (termStreamed) {

                lhs << (rhs.j() > 0 ? " + " : " - ");
                lhs << std::abs(rhs.j());

            } else {

                lhs << rhs.j();
                termStreamed = true;
            }

            lhs << "j";
        }

        if (!kZero) {

            if (termStreamed) {

                lhs << (rhs.k() > 0 ? " + " : " - ");
                lhs << std::abs(rhs.k());

            } else {

                lhs << rhs.k();
                termStreamed = true;
            }

            lhs << "k";
        }   

        if (multipleTerms) lhs << ")";

        return lhs;
    }

    // Alias types for coefficients of type int, long, float, double

    using iquat = tquat<int>;
    using lquat = tquat<long>;
    using fquat = tquat<float>;
    using dquat = tquat<double>;

    using quat = dquat;
}

namespace std {

    // Overload to call magn()
    template <typename T>
    constexpr double abs(const mth::tquat<T> &x) noexcept {

        return x.magn();
    }

    // Hash operator for use in certain STL containers
    template <typename T>
    struct hash<mth::tquat<T>> {

        size_t operator()(const mth::tquat<T> &x) {

            auto a = x.real();
            auto b = x.imaginary();

            return hash<decltype(a)>()(a) ^ hash<decltype(b)>()(b);
        }
    };
}

#endif
