#ifndef mth_tcomp_h__
#define mth_tcomp_h__

/* <mth/comp.h> - complex number header
 *      This includes the template class tcomp representing a complex number
 *      with coefficients of type T. Basic arithmetic operators are overloaded
 *      and member functions to find values such as the modulus and argument
 *      are defined. Functions for converting between cartesian and polar form
 *      are also provided along with more advanced math functions such as
 *      mth::exp and mth::cos.
 */

#include <cmath>
#include <ostream>
#include <iomanip>

#include <mth/mth.h>

namespace mth {

    // Forward declaration to avoid circular includes
    template <typename T, size_t N>
    class tvec;

    template <typename T>
    class tcomp {

    private:

        // Default initialize to zero
        T a = 0;
        T b = 0;

        // Direct initializations privatized to hide the implementation
        
        constexpr tcomp(const T &a, const T &b) noexcept
            :a(a), b(b) {}

        constexpr tcomp(const tvec<T, 2> &vec) noexcept
            :a(vec.x()), b(vec.y()) {}

    public:

        // Initialize to zero
        constexpr tcomp() noexcept = default;

        // Initialize to a real number
        constexpr tcomp(const T &a) noexcept
            :a(a) {}

        // Const and non-const getters
        // Non-const returns a reference for modification
#define BINDING(name, value) constexpr const T &name() const noexcept { return value; }\
                             constexpr T &name() noexcept { return value; }

        BINDING(real, a)
        BINDING(imag, b)

#undef BINDING

        // Casts distribute over the elements
        template <typename U>
        constexpr operator tcomp<U>() const noexcept {

            return tcomp<U>{static_cast<U>(a), static_cast<U>(b)};
        }

        // Convert to vector form
        constexpr tvec<T, 2> asCartesian() const noexcept {

            return tvec<T, 2>{a, b};
        }

        // Convert to principal polar form in a vector
        // Converts abs and arg from double to T for consistency
        constexpr tvec<T, 2> asPolar() const noexcept {

            return tvec<T, 2>{static_cast<T>(abs()), static_cast<T>(arg())};
        }

        // Returns the square magnitude, to avoid unnecessary calls to sqrt
        constexpr T absSqr() const noexcept {

            return a * a + b * b;
        }

        // Converted to double for more precision in sqrt 
        constexpr double abs() const noexcept {

            auto ls = static_cast<double>(absSqr());
            
            // If the magnitude is zero don't bother with sqrt
            if (util::isZero(ls)) return 0.0;

            return std::sqrt(ls);
        }

        // Returns the principal argument
        // Converted to double for more precision in atan2
        constexpr double arg() const noexcept {

            return std::atan2(static_cast<double>(b), static_cast<double>(a));
        }

        // Returns z / |z|
        constexpr tcomp<T> unit() const noexcept {

            return *this / abs();
        }

        // Returns the complex conjugate
        constexpr tcomp<T> conjugate() const noexcept {

            return tcomp<T>{a, -b};
        }

        // Returns 1 / z
        constexpr tcomp<T> inverse() const noexcept {

            return conjugate() / absSqr();
        }

        // Returns the complex representation of a rotation about the origin
        // Equivalent to fromPolar(1, angle)
        constexpr static tcomp<T> rotation(const T &angle) noexcept {

            // Convert to double for accurate cos/sin then convert back
            
            auto a = static_cast<double>(angle);

            auto c = static_cast<T>(std::cos(a));
            auto s = static_cast<T>(std::sin(a));

            return tcomp<T>{c, s};
        }

        // Create a complex number from its cartesian coordinates
        constexpr static tcomp<T> fromCartesian(const T& x, const T &y) noexcept {

            return tcomp<T>{x, y};
        }

        // Create a complex number from its cartesian vector
        constexpr static tcomp<T> fromCartesian(const tvec<T, 2> &vec) noexcept {

            return fromCartesian(vec.x(), vec.y());
        }

        // Create a complex number from its polar coordinates
        constexpr static tcomp<T> fromPolar(const T &radius, const T &angle) noexcept {

            return radius * rotation(angle);
        }

        // Create a complex number from its polar coordinates in a vector
        constexpr static tcomp<T> fromPolar(const tvec<T, 2> &polar) noexcept {

            return fromPolar(polar.x(), polar.y());
        }

        // Compound assignment operators

        constexpr tcomp<T> &operator+=(const tcomp<T> &rhs) noexcept {

            a += rhs.a;
            b += rhs.b;

            return *this;
        }

        constexpr tcomp<T> &operator+=(const T &rhs) noexcept {

            a += rhs;

            return *this;
        }

        constexpr tcomp<T> &operator-=(const tcomp<T> &rhs) noexcept {

            a -= rhs.a;
            b -= rhs.b;

            return *this;
        }

        constexpr tcomp<T> &operator-=(const T &rhs) noexcept {

            a -= rhs;

            return *this;
        }

        constexpr tcomp<T> &operator*=(const tcomp<T> &rhs) noexcept {

            *this = *this * rhs;

            return *this;
        }

        constexpr tcomp<T> &operator*=(const T &rhs) noexcept {

            a *= rhs;
            b *= rhs;

            return *this;
        }

        constexpr tcomp<T> &operator/=(const tcomp<T> &rhs) noexcept {

            *this = *this / rhs;

            return *this;
        }

        constexpr tcomp<T> &operator/=(const T &rhs) noexcept {

            a /= rhs;
            b /= rhs;

            return *this;
        }
    };

    // Alias types for coefficient types int, long, float, double 
    using icomp = tcomp<int>;
    using lcomp = tcomp<long>;
    using fcomp = tcomp<float>;
    using dcomp = tcomp<double>;
    using comp = dcomp;

    // The constant i for convenience
    template <typename T>
    constexpr tcomp<T> i = tcomp<T>::fromCartesian(0, 1);

    // Generic maths functions
    
    // Defers to tcomp<T>::abs
    template <typename T>
    constexpr double abs(const mth::tcomp<T> &z) noexcept {

        return z.abs();
    }

    // Calculate the principal sqrt of a complex number
    template <typename T>
    constexpr mth::tcomp<T> sqrt(const mth::tcomp<T> &z) noexcept {

        using std::sqrt;

        auto p = z.asPolar();

        p.x() = sqrt(p.x());
        p.y() /= 2;

        return mth::tcomp<T>::fromPolar(p);
    }

    // Calculate exp(z)
    template <typename T>
    constexpr mth::tcomp<T> exp(const mth::tcomp<T> &z) noexcept {

        using std::cos;
        using std::sin;
        using std::exp;

        auto c = cos(z.imag());
        auto s = sin(z.imag());

        return exp(z.real()) * (c + mth::i<T> * s);
    }

    // Calculate the principal log of a complex number
    template <typename T>
    constexpr mth::tcomp<T> log(const mth::tcomp<T> &z) noexcept {

        using std::log;

        auto p = z.asPolar();

        return log(p.x()) + mth::i<T> * p.y();
    }

    // Calculate cos(z)
    template <typename T>
    constexpr mth::tcomp<T> cos(const mth::tcomp<T> &z) noexcept {

        auto exponent = mth::i<T> * z;

        return (exp(exponent) + exp(-exponent)) / 2;
    }

    // Calculate sin(z)
    template <typename T>
    constexpr mth::tcomp<T> sin(const mth::tcomp<T> &z) noexcept {

        auto exponent = mth::i<T> * z;

        return (exp(exponent) - exp(-exponent)) / (2 * mth::i<T>);
    }

    // Calculate a complex power of a complex number using exp and log
    template <typename T>
    constexpr mth::tcomp<T> pow(const mth::tcomp<T> &z, const mth::tcomp<T> &exponent) noexcept {

        using std::log;

        return exp(exponent * log(z));
    }

    // Calculate a complex power of a real number using exp and log
    template <typename T>
    constexpr mth::tcomp<T> pow(const T &base, const mth::tcomp<T> &z) noexcept {

        using std::log;

        return exp(z * log(base));
    }

    // Calculate a real power of a complex number using exp and log
    template <typename T>
    constexpr mth::tcomp<T> pow(const mth::tcomp<T> &z, const T &exponent) noexcept {

        using std::log;

        return exp(exponent * log(z));
    }

    // Calculate a positive integer power of a complex number as a product
    template <typename T>
    constexpr mth::tcomp<T> pow(const mth::tcomp<T> &z, size_t exponent) noexcept {

        auto result = mth::tcomp<T>{1};

        for (size_t i = 1; i <= exponent; i++) {

            result *= z;
        }

        return result;
    }

    // TODO: Template for arbitrary scalars since implicit type conversion can't happen (and also in other classes)
    
    // Arithmetic operators

    template <typename T>
    constexpr tcomp<T> operator+(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result += rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator+(const tcomp<T> &lhs, const T &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result += rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator+(const T &lhs, const tcomp<T> &rhs) noexcept {

        return rhs + lhs;
    }

    template <typename T>
    constexpr tcomp<T> operator-(const tcomp<T> &rhs) noexcept {

        auto result = tcomp<T>{rhs};

        result.real() = -result.real();
        result.imag() = -result.imag();

        return result;
    }

    template <typename T>
    constexpr tcomp<T> operator-(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result -= rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator-(const tcomp<T> &lhs, const T &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result -= rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator-(const T &lhs, const tcomp<T> &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result -= rhs;
    }

    // Product is expanded as cartesian form
    template <typename T>
    constexpr tcomp<T> operator*(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        return tcomp<T>::fromCartesian(lhs.real() * rhs.real() - lhs.imag() * rhs.imag(),
                                       lhs.real() * rhs.imag() + lhs.imag() * rhs.real());
    }

    template <typename T>
    constexpr tcomp<T> operator*(const tcomp<T> &lhs, const T &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result *= rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator*(const T &lhs, const tcomp<T> &rhs) noexcept {

        return rhs * lhs;
    }

    template <typename T>
    constexpr tcomp<T> operator/(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        return lhs * rhs.inverse();
    }

    template <typename T>
    constexpr tcomp<T> operator/(const tcomp<T> &lhs, const T &rhs) noexcept {

        auto result = tcomp<T>{lhs};

        return result /= rhs;
    }

    template <typename T>
    constexpr tcomp<T> operator/(const T &lhs, const tcomp<T> &rhs) noexcept {

        return lhs * rhs.inverse();
    }

    // Equality operators use util::isEqual
    
    template <typename T>
    constexpr bool operator==(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        return util::isEqual(lhs.real(), rhs.real()) && util::isEqual(lhs.imag(), rhs.imag());
    }

    template <typename T>
    constexpr bool operator==(const tcomp<T> &lhs, const T &rhs) noexcept {

        return util::isZero(lhs.imag()) && util::isEqual(lhs.real(), rhs);
    }

    template <typename T>
    constexpr bool operator==(const T &lhs, const tcomp<T> &rhs) noexcept {

        return rhs == lhs;
    }

    template <typename T>
    constexpr bool operator!=(const tcomp<T> &lhs, const tcomp<T> &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T>
    constexpr bool operator!=(const tcomp<T> &lhs, const T &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T>
    constexpr bool operator!=(const T &lhs, const tcomp<T> &rhs) noexcept {

        return !(lhs == rhs);
    }

    // Pretty print
    template <typename T>
    std::ostream &operator<<(std::ostream &lhs, const tcomp<T> &rhs) {

        auto realZero = util::isZero(rhs.real());
        auto imagZero = util::isZero(rhs.imag());

        auto allZero = realZero && imagZero;
        auto noneZero = !realZero && !imagZero;

        if (allZero) return lhs << "0";

        // Bracket if there are multiple terms
        if (noneZero) lhs << "(";

        if (!realZero) lhs << rhs.real();

        if (!imagZero) {

            if (!realZero) {

                lhs << (rhs.imag() > 0 ? " + " : " - ");
                lhs << std::abs(rhs.imag());

            } else {

                lhs << rhs.imag();
            }

            lhs << "i";
        }

        if (noneZero) lhs << ")";

        return lhs;
    }
}

namespace std {

    // Hash function for use in certain STL containers
    template<typename T>
    struct hash<mth::tcomp<T>> {

        size_t operator()(const mth::tcomp<T> &z) const {

            auto r = z.real();
            auto i = z.imag();

            return hash<decltype(r)>()(r) ^ hash<decltype(i)>()(i);
        }
    };
}

#endif
