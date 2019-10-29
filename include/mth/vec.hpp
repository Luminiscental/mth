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
}
