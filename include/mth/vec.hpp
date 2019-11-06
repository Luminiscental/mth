#pragma once

#include <array>
#include <cstdint>

namespace mth
{
    // Vector with N elements of type T.
    template <typename T, size_t N>
    class tvec
    {
    private:
        // Components:

        std::array<T, N> _array;

        // Helper functions:

        template <typename Func, size_t... Ns>
        constexpr auto mapHelper(std::index_sequence<Ns...>, Func f) const
        {
            using O = decltype(f(get(0)));
            return tvec<O, N>{(f(get(Ns)))...};
        }

    public:
        // Constructors:

        constexpr tvec() : _array{} {}

        template <typename... Ts>
        constexpr tvec(Ts... elems) : _array{elems...}
        {
        }

        // Component accessors:

        constexpr T get(size_t index) const
        {
            return _array.at(index);
        }

        constexpr T operator[](size_t index) const
        {
            return get(index);
        }

        // Component aliases:

        constexpr T x() const
        {
            static_assert(N > 0, "vector does not have an x component");
            return get(0);
        }

        constexpr T y() const
        {
            static_assert(N > 1, "vector does not have a y component");
            return get(1);
        }

        constexpr T z() const
        {
            static_assert(N > 2, "vector does not have a z component");
            return get(2);
        }

        constexpr T w() const
        {
            static_assert(N > 3, "vector does not have a w component");
            return get(3);
        }

        constexpr T r() const
        {
            static_assert(N >= 3 && N <= 4, "vector is not a color size");
            return get(0);
        }

        constexpr T g() const
        {
            static_assert(N >= 3 && N <= 4, "vector is not a color size");
            return get(1);
        }

        constexpr T b() const
        {
            static_assert(N >= 3 && N <= 4, "vector is not a color size");
            return get(2);
        }

        constexpr T a() const
        {
            static_assert(N == 4, "vector is not a color with alpha size");
            return get(3);
        }

        // Iterators:

        constexpr auto begin() const noexcept
        {
            return _array.cbegin();
        }

        constexpr auto begin() noexcept
        {
            return _array.begin();
        }

        constexpr auto cbegin() const noexcept
        {
            return _array.cbegin();
        }

        constexpr auto end() noexcept
        {
            return _array.end();
        }

        constexpr auto end() const noexcept
        {
            return _array.cend();
        }

        constexpr auto cend() const noexcept
        {
            return _array.cend();
        }

        constexpr auto rbegin() noexcept
        {
            return _array.rbegin();
        }

        constexpr auto rbegin() const noexcept
        {
            return _array.rbegin();
        }

        constexpr auto crbegin() const noexcept
        {
            return _array.crbegin();
        }

        constexpr auto rend() noexcept
        {
            return _array.rend();
        }

        constexpr auto rend() const noexcept
        {
            return _array.rend();
        }

        constexpr auto crend() const noexcept
        {
            return _array.crend();
        }

        // Transformations / calculations:

        template <typename Func>
        constexpr auto map(Func f) const
        {
            return mapHelper(std::make_index_sequence<N>{}, f);
        }
    };

    // Arithmetic helper functions:

    template <typename T, size_t N, size_t... Ns>
    constexpr tvec<T, N>
    tvecSumHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return {(lhs.get(Ns) + rhs.get(Ns))...};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr tvec<T, N>
    tvecDiffHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return {(lhs.get(Ns) - rhs.get(Ns))...};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr tvec<T, N>
    tvecScaleHelper(std::index_sequence<Ns...>, T scalar, tvec<T, N> vector)
    {
        return {(scalar * vector.get(Ns))...};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr tvec<T, N>
    tvecReduceHelper(std::index_sequence<Ns...>, T scalar, tvec<T, N> vector)
    {
        return {(vector.get(Ns) / scalar)...};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr bool
    tvecEqualsHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return ((lhs.get(Ns) == rhs.get(Ns)) && ...);
    }

    // Operator overloads:

    template <typename T, size_t N>
    constexpr tvec<T, N> operator+(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return tvecSumHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator-(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return tvecDiffHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator*(T lhs, tvec<T, N> rhs)
    {
        return tvecScaleHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator*(tvec<T, N> lhs, T rhs)
    {
        return tvecScaleHelper(std::make_index_sequence<N>{}, rhs, lhs);
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator/(tvec<T, N> lhs, T rhs)
    {
        return tvecReduceHelper(std::make_index_sequence<N>{}, rhs, lhs);
    }

    template <typename T, size_t N>
    constexpr bool operator==(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return tvecEqualsHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }

    template <typename T, size_t N>
    constexpr bool operator!=(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return !(lhs == rhs);
    }

    // Aliases for N=2..4 with prefixes: - f for float
    //                                   - d for double
    //                                   - i for int
    //                                   - u for unsigned int

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

    // Aliases for default type float:

    using vec2 = fvec2;
    using vec3 = fvec3;
    using vec4 = fvec4;
}
