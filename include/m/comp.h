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

    template <typename T>
    auto operator+(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator+(const tcomp<T> &lhs, const T &rhs); 

    template <typename T>
    auto operator+(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator-(const tcomp<T> &rhs);

    template <typename T>
    auto operator-(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator-(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    auto operator-(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator*(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator*(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    auto operator*(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator/(const tcomp<T> &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto operator/(const tcomp<T> &lhs, const T &rhs);

    template <typename T>
    auto operator/(const T &lhs, const tcomp<T> &rhs);

    template <typename T>
    auto &operator<<(std::ostream &lhs, const tcomp<T> &rhs);

    template <typename T>
    class tcomp {

    private:

        T a;
        T b;

    public:

        constexpr tcomp() noexcept
            :a(0), b(0) {}

        constexpr tcomp(const T &a) noexcept
            :a(a), b(0) {}

        constexpr tcomp(const T &a, const T &b) noexcept
            :a(a), b(b) {}

        constexpr tcomp(tvec<T, 2> &vec) noexcept
            :a(vec.x()), b(vec.y()) {}

#define BINDING(name, value)    const auto & name () const noexcept; \
                                auto & name () noexcept;

        BINDING(real, a)
        BINDING(imag, b)

#undef BINDING

        tvec<T, 2> asPolar() const noexcept;

        operator tvec<T, 2>() noexcept;

        auto absSqr() const noexcept;

        auto abs() const noexcept;

        auto arg() const noexcept;

        auto unit() const;

        auto conjugate() const noexcept;

        auto inverse() const;

        static tcomp<T> rotation(const T &angle);

        static tcomp<T> fromPolar(const T &radius, const T &angle);

        static tcomp<T> fromPolar(const tvec<T, 2> &polar);

        auto &operator+=(const tcomp<T> &rhs);

        auto &operator+=(const T &rhs);

        auto &operator-=(const tcomp<T> &rhs);

        auto &operator-=(const T &rhs);

        auto &operator*=(const tcomp <T> &rhs);

        auto &operator*=(const T &rhs);

        auto &operator/=(const tcomp<T> &rhs);

        auto &operator/=(const T &rhs);
    };

    using icomp = tcomp<int>;
    using lcomp = tcomp<long>;
    using  comp = tcomp<float>;
    using dcomp = tcomp<double>;

    // TODO: Should this be a template?
    constexpr tcomp<float> i(0, 1);

    auto operator "" _i(long double im);
}

namespace std {

    template <typename T>
    double abs(const m::tcomp<T> &z);

    template <typename T>
    m::tcomp<T> sqrt(const m::tcomp<T> &z);

    template <typename T>
    m::tcomp<T> exp(const m::tcomp<T> &z);

    template <typename T>
    m::tcomp<T> cos(const m::tcomp<T> &z);

    template <typename T>
    m::tcomp<T> sin(const m::tcomp<T> &z);
}

#include <m/comp_impl.h>

#endif
