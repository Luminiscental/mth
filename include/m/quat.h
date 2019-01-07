#ifndef __m_quat_h__
#define __m_quat_h__

#include <ostream>

#include <m/vec.h>

namespace m {

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

        T magnSqr() const noexcept;

        double magn() const noexcept;

        constexpr tquat<T> conjugate() const noexcept {

            return tquat<T>(w, -ijk);
        }

        tquat<T> inverse() const;

        tquat<T> unit() const;

        tvec<T, 3> rotate(const m::tvec<T, 3> &vector) const;

        static tquat<T> identity();

                                                               //                        /-
        static tquat<T> rotation(T angle, const tvec<T, 3> &axis); // NOTE: Right-handed: ---|--> axis
                                                               //                        \-> rotation

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

    using iquat = tquat<int>;
    using lquat = tquat<long>;
    using  quat = tquat<float>;
    using dquat = tquat<double>;
}

#include <m/quat_impl.h>

#endif
