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
        constexpr tcomp()
        {
        }

        template <typename Comp, typename = enable_if_comp_t<Comp>>
        constexpr tcomp(Comp expr) : _real{expr.real()}, _imag{expr.imag()}
        {
        }

        constexpr tcomp(T real, T imag) : _real{real}, _imag{imag} {}

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
    };

    using fcomp = tcomp<float>;
    using dcomp = tcomp<double>;
    using icomp = tcomp<int>;
    using ucomp = tcomp<unsigned int>;

    using comp = fcomp;
}
