#pragma once

namespace mth
{
    namespace util
    {
        // Given an arithmetic type T, tries to give the representation of zero
        // in that type.
        template <typename T>
        [[nodiscard]] constexpr T zero() noexcept
        {
            return T{} * 0;
        }

        // Given a value of an arithemtic type T, checks for equality with
        // mth::util::zero<T>()
        template <typename T>
        [[nodiscard]] constexpr bool isZero(T value) noexcept
        {
            return value == zero<T>();
        }

        // Given a value and min/max of an arithemtic type T, checks if the
        // value is within the range including equality at the endpoints.
        template <typename T>
        [[nodiscard]] constexpr bool inRangeClosed(T value, T min, T max)
        {
            return min <= value && value <= max;
        }

        // Given a value and min/max of an arithemtic type T, checks if the
        // value is within the range excluding equality at the endpoints.
        template <typename T>
        [[nodiscard]] constexpr bool inRangeOpen(T value, T min, T max)
        {
            return min < value && value < max;
        }
    }
}
