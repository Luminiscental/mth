#pragma once

#include <array>
#include <cstdint>

namespace mth
{
    template <typename Derived>
    class tcomp_expr;

    // tcomp_info declarations

    template <typename Comp>
    class tcomp_info;

    template <typename Comp>
    using tcomp_elem_t = typename tcomp_info<Comp>::Elem;

    template <typename... Ts>
    using enable_if_comps_t =
        std::enable_if_t<((std::is_base_of_v<tcomp_expr<Ts>, Ts>) &&...)>;

    template <typename T>
    using enable_if_comp_t = enable_if_comps_t<T>;

    // tcomp definition

    template <typename T>
    class tcomp : public tcomp_expr<tcomp<T>>
    {
    private:
        T _real;
        T _imag;

    public:
        template <
            typename = std::enable_if_t<std::is_default_constructible_v<T>>>
        constexpr tcomp() : _real{}, _imag{}
        {
        }

        template <typename Comp, typename = enable_if_comp_t<Comp>>
        constexpr tcomp(Comp expr) : _real{expr.real()}, _imag{expr.imag()}
        {
        }

        constexpr tcomp(T real, T imag)
            : _real{std::move(real)}, _imag{std::move(imag)}
        {
        }

        constexpr T real() const
        {
            return _real;
        }

        constexpr T imag() const
        {
            return _imag;
        }
    };

    template <typename T>
    class tcomp_info<tcomp<T>>
    {
    public:
        using Elem = T;
    };

    // tcomp_sum definition

    template <typename Lhs, typename Rhs>
    class tcomp_sum : public tcomp_expr<tcomp_sum<Lhs, Rhs>>
    {
    private:
        Lhs _lhs;
        Rhs _rhs;

    public:
        explicit constexpr tcomp_sum(Lhs lhs, Rhs rhs)
            : _lhs{std::move(lhs)}, _rhs{std::move(rhs)}
        {
        }

        constexpr tcomp_elem_t<Lhs> real() const
        {
            return _lhs.real() + _rhs.real();
        }

        constexpr tcomp_elem_t<Lhs> imag() const
        {
            return _lhs.imag() + _rhs.imag();
        }
    };

    template <typename Lhs, typename Rhs>
    class tcomp_info<tcomp_sum<Lhs, Rhs>>
    {
    public:
        using Elem = tcomp_elem_t<Lhs>;
    };

    // tcomp_neg definition

    template <typename Comp>
    class tcomp_neg : public tcomp_expr<tcomp_neg<Comp>>
    {
    private:
        Comp _target;

    public:
        explicit constexpr tcomp_neg(Comp comp) : _target{std::move(comp)} {}

        constexpr tcomp_elem_t<Comp> real() const
        {
            return -_target.real();
        }

        constexpr tcomp_elem_t<Comp> imag() const
        {
            return -_target.imag();
        }
    };

    template <typename Comp>
    class tcomp_info<tcomp_neg<Comp>>
    {
    public:
        using Elem = tcomp_elem_t<Comp>;
    };

    // tcomp_conj definition

    template <typename Comp>
    class tcomp_conj : public tcomp_expr<tcomp_conj<Comp>>
    {
    private:
        Comp _target;

    public:
        explicit constexpr tcomp_conj(Comp comp) : _target{std::move(comp)} {}

        constexpr tcomp_elem_t<Comp> real() const
        {
            return _target.real();
        }

        constexpr tcomp_elem_t<Comp> imag() const
        {
            return -_target.imag();
        }
    };

    template <typename Comp>
    class tcomp_info<tcomp_conj<Comp>>
    {
    public:
        using Elem = tcomp_elem_t<Comp>;
    };

    // tcomp_expr definition

    template <typename Derived>
    class tcomp_expr
    {
    public:
        using Elem = tcomp_elem_t<Derived>;

    protected:
        // CRTP safeguard
        tcomp_expr()  = default;
        ~tcomp_expr() = default;

    public:
        // Component accessors:

        constexpr Elem real() const
        {
            return static_cast<Derived const&>(*this).real();
        }

        constexpr Elem imag() const
        {
            return static_cast<Derived const&>(*this).imag();
        }

        constexpr auto conj() const
        {
            Derived copy = static_cast<Derived const&>(*this);
            return tcomp_conj{std::move(copy)};
        }
    };

    enum class tcomp_tag
    {
        Dummy
    };

    // operator overloads

    template <
        typename Lhs,
        typename Rhs,
        typename  = enable_if_comps_t<Lhs, Rhs>,
        tcomp_tag = tcomp_tag::Dummy>
    constexpr auto operator+(Lhs lhs, Rhs rhs)
    {
        return tcomp_sum{std::move(lhs), std::move(rhs)};
    }

    template <
        typename Comp,
        typename  = enable_if_comp_t<Comp>,
        tcomp_tag = tcomp_tag::Dummy>
    constexpr auto operator-(Comp comp)
    {
        return tcomp_neg{std::move(comp)};
    }

    // aliases

    using fcomp = tcomp<float>;
    using dcomp = tcomp<double>;
    using icomp = tcomp<int>;
    using ucomp = tcomp<unsigned int>;

    using comp = fcomp;
}
