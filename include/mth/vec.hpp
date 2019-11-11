#pragma once

#include <array>
#include <cstdint>

namespace mth
{
    template <typename Derived>
    class tvec_expr;

    // tvec_info declarations

    template <typename Vec>
    class tvec_info;

    template <typename Vec>
    using elem_t = typename tvec_info<Vec>::Elem;

    template <typename Vec>
    constexpr auto size_v = tvec_info<Vec>::SIZE;

    template <typename... Ts>
    using enable_if_vecs_t =
        std::enable_if_t<((std::is_base_of_v<tvec_expr<Ts>, Ts>) &&...)>;

    template <typename T>
    using enable_if_vec_t = enable_if_vecs_t<T>;

    // tvec definition

    template <typename T, size_t N>
    class tvec : public tvec_expr<tvec<T, N>>
    {
    private:
        std::array<T, N> _components;

        template <typename Vec, size_t... Ns>
        tvec(std::index_sequence<Ns...>, Vec expr)
            : _components{expr.get(Ns)...}
        {
        }

    public:
        constexpr tvec() : _components{} {}

        template <typename Vec, typename = enable_if_vec_t<Vec>>
        tvec(Vec expr)
            : tvec{std::make_index_sequence<size_v<Vec>>{}, std::move(expr)}
        {
            static_assert(
                size_v<Vec> == N, "incompatible vector expression size");
        }

        template <typename... Ts>
        constexpr tvec(Ts... values) : _components{std::move(values)...}
        {
            static_assert(
                sizeof...(Ts) == N,
                "wrong number of components to initialize vector");
        }

        constexpr T get(size_t index) const
        {
            return _components.at(index);
        }

        // iterators (here rather than in tvec_expr since they require the
        // vector to be concrete)

        constexpr auto begin()
        {
            return _components.begin();
        }

        constexpr auto begin() const
        {
            return _components.begin();
        }

        constexpr auto cbegin() const
        {
            return _components.cbegin();
        }

        constexpr auto rbegin()
        {
            return _components.rbegin();
        }

        constexpr auto crbegin() const
        {
            return _components.crbegin();
        }

        constexpr auto end()
        {
            return _components.end();
        }

        constexpr auto end() const
        {
            return _components.end();
        }

        constexpr auto cend() const
        {
            return _components.cend();
        }

        constexpr auto rend()
        {
            return _components.rend();
        }

        constexpr auto crend() const
        {
            return _components.crend();
        }
    };

    template <typename Vec, typename = enable_if_vec_t<Vec>>
    tvec(Vec expr)->tvec<elem_t<Vec>, size_v<Vec>>;

    template <typename T, size_t N>
    class tvec_info<tvec<T, N>>
    {
    public:
        using Elem                 = T;
        static constexpr auto SIZE = N;
    };

    // tvec_sum definition

    template <typename Lhs, typename Rhs>
    class tvec_sum : public tvec_expr<tvec_sum<Lhs, Rhs>>
    {
    private:
        Lhs _lhs;
        Rhs _rhs;

    public:
        explicit constexpr tvec_sum(Lhs lhs, Rhs rhs)
            : _lhs{std::move(lhs)}, _rhs{std::move(rhs)}
        {
        }

        constexpr elem_t<Lhs> get(size_t index) const
        {
            return _lhs.get(index) + _rhs.get(index);
        }
    };

    template <typename Lhs, typename Rhs>
    class tvec_info<tvec_sum<Lhs, Rhs>>
    {
    public:
        using Elem                 = elem_t<Lhs>;
        static constexpr auto SIZE = size_v<Lhs>;
    };

    // tvec_diff definition

    template <typename Lhs, typename Rhs>
    class tvec_diff : public tvec_expr<tvec_diff<Lhs, Rhs>>
    {
    private:
        Lhs _lhs;
        Rhs _rhs;

    public:
        explicit constexpr tvec_diff(Lhs lhs, Rhs rhs)
            : _lhs{std::move(lhs)}, _rhs{std::move(rhs)}
        {
        }

        constexpr elem_t<Lhs> get(size_t index) const
        {
            return _lhs.get(index) - _rhs.get(index);
        }
    };

    template <typename Lhs, typename Rhs>
    class tvec_info<tvec_diff<Lhs, Rhs>>
    {
    public:
        using Elem                 = elem_t<Lhs>;
        static constexpr auto SIZE = size_v<Lhs>;
    };

    // tvec_scale definition

    template <typename T, typename Vec>
    class tvec_scale : public tvec_expr<tvec_scale<T, Vec>>
    {
    private:
        Vec _vector;
        T _scalar;

    public:
        explicit constexpr tvec_scale(T scalar, Vec vector)
            : _vector{std::move(vector)}, _scalar{std::move(scalar)}
        {
        }

        constexpr T get(size_t index) const
        {
            return _scalar * _vector.get(index);
        }
    };

    template <typename T, typename Vec>
    class tvec_info<tvec_scale<T, Vec>>
    {
    public:
        using Elem                 = elem_t<Vec>;
        static constexpr auto SIZE = size_v<Vec>;
    };

    // tvec_reduce definition

    template <typename T, typename Vec>
    class tvec_reduce : public tvec_expr<tvec_reduce<T, Vec>>
    {
    private:
        Vec _vector;
        T _scalar;

    public:
        explicit constexpr tvec_reduce(T scalar, Vec vector)
            : _vector{std::move(vector)}, _scalar{std::move(scalar)}
        {
        }

        constexpr T get(size_t index) const
        {
            return _vector.get(index) / _scalar;
        }
    };

    template <typename T, typename Vec>
    class tvec_info<tvec_reduce<T, Vec>>
    {
    public:
        using Elem                 = elem_t<Vec>;
        static constexpr auto SIZE = size_v<Vec>;
    };

    // note: Vecs must be non-empty as std::tuple requires so

    template <typename Func, typename... Vecs>
    class tvec_map : public tvec_expr<tvec_map<Func, Vecs...>>
    {
    private:
        std::tuple<Vecs...> _vecs;
        Func _functor;

        template <size_t... Ns>
        constexpr auto getHelper(std::index_sequence<Ns...>, size_t index) const
        {
            return _functor((std::get<Ns>(_vecs).get(index))...);
        }

    public:
        explicit constexpr tvec_map(Func functor, Vecs... vecs)
            : _vecs{vecs...}, _functor{std::move(functor)}
        {
        }

        constexpr auto get(size_t index) const
        {
            return getHelper(
                std::make_index_sequence<sizeof...(Vecs)>{}, index);
        }
    };

    // note: only specialized case defined since empty parameter pack isn't
    // valid as above

    template <typename Func, typename Vec, typename... Vecs>
    class tvec_info<tvec_map<Func, Vec, Vecs...>>
    {
    public:
        using Elem = std::result_of_t<Func(elem_t<Vec>, elem_t<Vecs>...)>;
        static constexpr auto SIZE = size_v<Vec>;
    };

    // tvec_expr definition

    // TODO: memoize?

    template <typename Derived>
    class tvec_expr
    {
    protected:
        // CRTP safeguard
        tvec_expr()  = default;
        ~tvec_expr() = default;

    public:
        using Elem                 = elem_t<Derived>;
        static constexpr auto SIZE = size_v<Derived>;

        // Component accessors:

        constexpr Elem get(size_t index) const
        {
            return static_cast<Derived const*>(this)->get(index);
        }

        constexpr Elem operator[](size_t index) const
        {
            return get(index);
        }

        // component aliases:

        constexpr Elem x() const
        {
            static_assert(SIZE > 0, "vector does not have an x component");
            return get(0);
        }

        constexpr Elem y() const
        {
            static_assert(SIZE > 1, "vector does not have a y component");
            return get(1);
        }

        constexpr Elem z() const
        {
            static_assert(SIZE > 2, "vector does not have a z component");
            return get(2);
        }

        constexpr Elem w() const
        {
            static_assert(SIZE > 3, "vector does not have a w component");
            return get(3);
        }

        constexpr Elem r() const
        {
            static_assert(SIZE >= 3 && SIZE <= 4, "vector is not a color size");
            return get(0);
        }

        constexpr Elem g() const
        {
            static_assert(SIZE >= 3 && SIZE <= 4, "vector is not a color size");
            return get(1);
        }

        constexpr Elem b() const
        {
            static_assert(SIZE >= 3 && SIZE <= 4, "vector is not a color size");
            return get(2);
        }

        constexpr Elem a() const
        {
            static_assert(SIZE == 4, "vector is not a color with alpha size");
            return get(3);
        }

        // transforms / calculations

        template <typename Func>
        constexpr auto map(Func f) const
        {
            return tvec_map{f, static_cast<Derived const&>(*this)};
        }
    };

    // static functions

    namespace vec
    {
        template <
            typename Func,
            typename... Vecs,
            typename = enable_if_vecs_t<Vecs...>>
        constexpr auto map(Func f, Vecs... vecs)
        {
            return tvec_map{f, vecs...};
        }
    }

    template <typename Derived>
    class tvec_info<tvec_expr<Derived>>
    {
    public:
        using Elem                 = elem_t<Derived>;
        static constexpr auto SIZE = size_v<Derived>;
    };

    // operator overloads

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vec_t<Lhs>,
        typename = enable_if_vec_t<Rhs>>
    constexpr auto operator+(Lhs lhs, Rhs rhs)
    {
        return tvec_sum{lhs, rhs};
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vec_t<Lhs>,
        typename = enable_if_vec_t<Rhs>>
    constexpr auto operator-(Lhs lhs, Rhs rhs)
    {
        return tvec_diff{lhs, rhs};
    }

    template <typename T, typename Vec, typename = enable_if_vec_t<Vec>>
    constexpr auto operator*(T lhs, Vec rhs)
    {
        return tvec_scale{lhs, rhs};
    }

    template <typename T, typename Vec, typename = enable_if_vec_t<Vec>>
    constexpr auto operator*(Vec lhs, T rhs)
    {
        return rhs * lhs;
    }

    template <typename T, typename Vec, typename = enable_if_vec_t<Vec>>
    constexpr auto operator/(Vec lhs, T rhs)
    {
        return tvec_reduce{rhs, lhs};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr auto
    equalityHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return ((lhs.get(Ns) == rhs.get(Ns)) && ...);
    }

    template <typename T, size_t N>
    constexpr auto operator==(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return equalityHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }

    template <typename T, size_t N>
    constexpr auto operator!=(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return !(lhs == rhs);
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vec_t<Lhs>,
        typename = enable_if_vec_t<Rhs>>
    constexpr auto operator==(Lhs lhs, Rhs rhs)
    {
        return tvec{lhs} == tvec{rhs};
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vec_t<Lhs>,
        typename = enable_if_vec_t<Rhs>>
    constexpr auto operator!=(Lhs lhs, Rhs rhs)
    {
        return !(lhs == rhs);
    }

    // aliases

    using fvec2 = tvec<float, 2>;
    using fvec3 = tvec<float, 3>;
    using fvec4 = tvec<float, 4>;

    using dvec2 = tvec<double, 2>;
    using dvec3 = tvec<double, 3>;
    using dvec4 = tvec<double, 4>;

    using ivec2 = tvec<int, 2>;
    using ivec3 = tvec<int, 3>;
    using ivec4 = tvec<int, 4>;

    using uvec2 = tvec<unsigned int, 2>;
    using uvec3 = tvec<unsigned int, 3>;
    using uvec4 = tvec<unsigned int, 4>;

    using vec2 = fvec2;
    using vec3 = fvec3;
    using vec4 = fvec4;
}
