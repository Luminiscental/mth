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
        tvec() : _array{} {}

        template <typename... Ts>
        tvec(Ts... elems) : _array{elems...}
        {
        }

        T get(size_t index)
        {
            return _array.at(index);
        }
    };

    template <typename T, size_t N, size_t... Ns>
    tvec<T, N>
    tvecSumHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return tvec<T, N>{(lhs.get(Ns) + rhs.get(Ns))...};
    }

    template <typename T, size_t N>
    constexpr tvec<T, N> operator+(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return tvecSumHelper(std::make_index_sequence<N>{}, lhs, rhs);
    }
}
