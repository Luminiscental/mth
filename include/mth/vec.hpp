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
        std::array<T, N> _array;

    public:
        constexpr tvec() : _array{} {}

        template <typename... Ts>
        constexpr tvec(Ts... elems) : _array{elems...}
        {
        }

        constexpr T get(size_t index)
        {
            return _array.at(index);
        }

        constexpr T operator[](size_t index)
        {
            return get(index);
        }

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
    };

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

    // Aliases for default type is float

    using vec2 = fvec2;
    using vec3 = fvec3;
    using vec4 = fvec4;
}
