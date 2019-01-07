#ifndef __m_tcomp_h__
#define __m_tcomp_h__

#include <cmath>
#include <ostream>
#include <iomanip>

namespace m {

    // NOTE: Forward declaration to avoid circular includes
    template <typename T, size_t N>
    class tvec; 

    template <typename T>
    class tcomp;

    // TODO: Template for arbitrary scalars since implicit type conversion can't happen

    template <typename T>
    tcomp<T> operator+(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator+(const tcomp<T> &lhs, const T &rhs); 

    template <typename T>
    tcomp<T> operator+(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator-(const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator-(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator-(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    tcomp<T> operator-(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator*(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator*(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    tcomp<T> operator*(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator/(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    tcomp<T> operator/(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    tcomp<T> operator/(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    bool operator==(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    bool operator==(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    bool operator==(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    bool operator!=(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    bool operator!=(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    bool operator!=(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    std::ostream &operator<<(std::ostream &lhs, const tcomp<T> &rhs);

    template <typename T>
    class tcomp {

    private:

        T a;
        T b;

        constexpr tcomp(const T &a, const T &b) noexcept
            :a(a), b(b) {}

        constexpr tcomp(tvec<T, 2> &vec) noexcept
            :a(vec.x()), b(vec.y()) {}

    public:

        constexpr tcomp() noexcept
            :a(0), b(0) {}

        constexpr tcomp(const T &a) noexcept
            :a(a), b(0) {}

#define BINDING(name, value)    const T & name () const noexcept; \
                                T & name () noexcept;

        BINDING(real, a)
        BINDING(imag, b)

#undef BINDING

        tvec<T, 2> asCartesian() const noexcept;

        tvec<T, 2> asPolar() const noexcept;

        operator tvec<T, 2>() noexcept;

        T absSqr() const noexcept;

        double abs() const noexcept;

        double arg() const noexcept;

        tcomp<T> unit() const;

        tcomp<T> conjugate() const noexcept;

        tcomp<T> inverse() const;

        static tcomp<T> rotation(const T &angle);

        constexpr static tcomp<T> fromCartesian(const T &x, const T &y) {

            return tcomp<T>(x, y);
        }

        static tcomp<T> fromCartesian(const tvec<T, 2> &vec);

        static tcomp<T> fromPolar(const T &radius, const T &angle);

        static tcomp<T> fromPolar(const tvec<T, 2> &polar);

        tcomp<T> &operator+=(const tcomp<T> &rhs);

        tcomp<T> &operator+=(const T &rhs);

        tcomp<T> &operator-=(const tcomp<T> &rhs);

        tcomp<T> &operator-=(const T &rhs);

        tcomp<T> &operator*=(const tcomp <T> &rhs);

        tcomp<T> &operator*=(const T &rhs);

        tcomp<T> &operator/=(const tcomp<T> &rhs);

        tcomp<T> &operator/=(const T &rhs);

        template <typename U>
        friend tcomp<U> operator+(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator+(const tcomp<U> &lhs, const U &rhs); 

        template <typename U>
        friend tcomp<U> operator+(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator-(const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator-(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator-(const tcomp<U> &lhs, const U &rhs);

        template <typename U>
        friend tcomp<U> operator-(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator*(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator*(const tcomp<U> &lhs, const U &rhs);

        template <typename U>
        friend tcomp<U> operator*(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator/(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend tcomp<U> operator/(const tcomp<U> &lhs, const U &rhs);

        template <typename U>
        friend tcomp<U> operator/(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend bool operator==(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend bool operator==(const tcomp<U> &lhs, const U &rhs);

        template <typename U>
        friend bool operator==(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend bool operator!=(const tcomp<U> &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend bool operator!=(const tcomp<U> &lhs, const U &rhs);

        template <typename U>
        friend bool operator!=(const U &lhs, const tcomp<U> &rhs);

        template <typename U>
        friend std::ostream &operator<<(std::ostream &lhs, const tcomp<U> &rhs);
    };

    using icomp = tcomp<int>;
    using lcomp = tcomp<long>;
    using  comp = tcomp<float>;
    using dcomp = tcomp<double>;

    template <typename T>
    constexpr tcomp<T> i = tcomp<T>::fromCartesian(0, 1);
}

namespace std {

    // TODO: Figure out how to template these while still getting overload detection

    double abs(const m::comp &z);

    m::comp sqrt(const m::comp &z);

    m::comp exp(const m::comp &z);

    m::comp cos(const m::comp &z);

    m::comp sin(const m::comp &z);
    
    // TODO: Add for other types
    template<>
    struct hash<m::comp> {

        size_t operator()(const m::comp &x) const;
    };
}

#include <m/comp_impl.h>

#endif
