#ifndef __m_tcomp_h__
#define __m_tcomp_h__

/* <m/comp.h> - complex number header
 *      This includes the template class tcomp representing a complex number with coefficients of type T. Basic arithmetic
 *      operators are overloaded and member functions to find values such as the modulus and argument are defined and also
 *      functions for converting between cartesian and polar form. More advanced math functions such as std::exp and std::cos
 *      are also overloaded.
 */

#include <cmath>
#include <ostream>
#include <iomanip>

namespace m {

    // Forward declaration to avoid circular includes
    template <typename T, size_t N>
    class tvec; 

    // Forward declaration for friending

    template <typename T>
    class tcomp;

    // TODO: Template for arbitrary scalars since implicit type conversion can't happen (and also in other classes)

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

        // Direct initialization privatized to hide the implementation

        constexpr tcomp(const T &a, const T &b) noexcept
            :a(a), b(b) {}

        constexpr tcomp(tvec<T, 2> &vec) noexcept
            :a(vec.x()), b(vec.y()) {}

    public:

        // Default initializes to zero
        constexpr tcomp() noexcept
            :a(0), b(0) {}

        // Initialize to a real number
        constexpr tcomp(const T &a) noexcept
            :a(a), b(0) {}

#define BINDING(name, value)    const T & name () const noexcept; \
                                      T & name ()       noexcept;

        BINDING(real, a)
        BINDING(imag, b)

#undef BINDING

        // Convert to vector forms
        tvec<T, 2> asCartesian() const noexcept;

        // Converts abs and arg from double to T so possible loss of information
        tvec<T, 2> asPolar() const noexcept;

        T absSqr() const noexcept;

        // Converted to double for more precision in sqrt / atan2
        double abs() const noexcept;
        double arg() const noexcept;

        // Returns z / |z|
        tcomp<T> unit() const;

        tcomp<T> conjugate() const noexcept;

        // Returns 1 / z
        tcomp<T> inverse() const;

        // Returns the complex representation of a rotation about the origin in 2-dimensions (CCW)
        // - equivalent to fromPolar(1, angle)
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

    // Alias types for coefficient types int, long, float, double 

    using icomp = tcomp<int>;
    using lcomp = tcomp<long>;
    using  comp = tcomp<float>;
    using dcomp = tcomp<double>;

    // The constant i for convenience

    template <typename T>
    constexpr tcomp<T> i = tcomp<T>::fromCartesian(0, 1);
}

namespace std {

    // Overloads of common math functions in <cmath>

    // TODO: Figure out how to template these while still getting overload detection

    double abs(const m::comp &z);

    m::comp sqrt(const m::comp &z);

    m::comp exp(const m::comp &z);

    m::comp cos(const m::comp &z);

    m::comp sin(const m::comp &z);

    // Hash operator for use in certain STL containers
    
    // TODO: Add for other types
    template<>
    struct hash<m::comp> {

        size_t operator()(const m::comp &x) const;
    };
}

#include <m/comp_impl.h>

#endif
