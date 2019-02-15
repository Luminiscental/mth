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

    // Forward declaration for friending

    template <typename T>
    class tquat;

    template <typename T>
    tquat<T> operator+(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator+(const tquat<T> &lhs, T rhs);

    template <typename T>
    tquat<T> operator+(T lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator-(const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator-(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator-(const tquat<T> &lhs, T rhs);

    template <typename T>
    tquat<T> operator-(T lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator*(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator*(T lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator*(const tquat<T> &lhs, T rhs);

    template <typename T>
    tquat<T> operator/(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator/(T lhs, const tquat<T> &rhs);

    template <typename T>
    tquat<T> operator/(const tquat<T> &lhs, T rhs);

    template <typename T>
    bool operator==(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    bool operator==(const T &lhs, const tquat<T> &rhs);

    template <typename T>
    bool operator==(const tquat<T> &lhs, const T &rhs);

    template <typename T>
    bool operator!=(const tquat<T> &lhs, const tquat<T> &rhs);

    template <typename T>
    bool operator!=(const T &lhs, const tquat<T> &rhs);

    template <typename T>
    bool operator!=(const tquat<T> &lhs, const T &rhs);

    template <typename T>
    std::ostream &operator<<(std::ostream &lhs, const tquat<T> &rhs);

    template <typename T>
    class tquat {

    private:

        T w;
        tvec<T, 3> ijk;

    public:

#define BINDING(name, value)    const auto & name () const; \
                                      auto & name ();

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

        T magnSqr() const noexcept;

        // Converted to double for more accurate sqrt
        double magn() const noexcept;

        constexpr tquat<T> conjugate() const noexcept {

            return tquat<T>(w, -ijk);
        }

        // Returns the reciprocal
        tquat<T> inverse() const;

        tquat<T> unit() const;

        // Returns the rotation represented by *this applied to the given vector
        tvec<T, 3> rotate(const mth::tvec<T, 3> &vector) const;

        // Returns 1
        static tquat<T> identity();

        // Converts euler angle and axis to quaternion representation
        static tquat<T> rotation(T angle, const tvec<T, 3> &axis);

        tquat<T> &operator+=(const tquat<T> &rhs);
        tquat<T> &operator+=(T rhs);
        tquat<T> &operator-=(const tquat<T> &rhs);
        tquat<T> &operator-=(T rhs);
        tquat<T> &operator*=(const tquat<T> &rhs);
        tquat<T> &operator*=(T rhs);
        tquat<T> &operator/=(const tquat<T> &rhs);
        tquat<T> &operator/=(T rhs);

        template <typename U>
        friend tquat<U> operator+(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator+(const tquat<U> &lhs, U rhs);

        template <typename U>
        friend tquat<U> operator+(U lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator-(const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator-(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator-(const tquat<U> &lhs, U rhs);

        template <typename U>
        friend tquat<U> operator-(U lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator*(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator*(U lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator*(const tquat<U> &lhs, U rhs);

        template <typename U>
        friend tquat<U> operator/(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator/(U lhs, const tquat<U> &rhs);

        template <typename U>
        friend tquat<U> operator/(const tquat<U> &lhs, U rhs);

        template <typename U>
        friend bool operator==(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend bool operator==(const U &lhs, const tquat<U> &rhs);

        template <typename U>
        friend bool operator==(const tquat<U> &lhs, const U &rhs);

        template <typename U>
        friend bool operator!=(const tquat<U> &lhs, const tquat<U> &rhs);

        template <typename U>
        friend bool operator!=(const U &lhs, const tquat<U> &rhs);

        template <typename U>
        friend bool operator!=(const tquat<U> &lhs, const U &rhs);

        template <typename U>
        friend std::ostream &operator<<(std::ostream &lhs, const tquat<U> &rhs);
    };

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
    double abs(const mth::tquat<T> &x);

    // Hash operator for use in certain STL containers
    template <typename T>
    struct hash<mth::tquat<T>> {

        size_t operator()(const mth::tquat<T> &x);
    };
}

// TODO!: Refactor to remove this
#include <mth/quat_impl.h>

#endif
