#pragma once

/// @file mth/mth.hpp Header containing definitions for the `mth::util`
/// namespace.

namespace mth
{
    /** @brief Contains general mathematical utility functions. */
    namespace util
    {
        /**
         * @brief Given an arithmetic type `T`, tries to give the representation
         * of zero in that type.
         *
         * @tparam T The type to get 0 for.
         *
         * @return The literal 0 statically casted to `T`.
         */
        template <typename T>
        [[nodiscard]] constexpr T zero() noexcept
        {
            return static_cast<T>(0);
        }

        /**
         * @brief Given a value of an arithemtic type `T`, checks for equality
         * with `mth::util::zero<T>()`.
         *
         * @tparam T The value type to check for zero.
         * @param value The value to check for equality with zero.
         * @return The result of comparing for equality against
         * `mth::zero<T>()`.
         */
        template <typename T>
        [[nodiscard]] constexpr bool isZero(T value) noexcept
        {
            return value == zero<T>();
        }

        /**
         * @brief Given a value and min/max of an arithemtic type `T`, checks if
         * the value is within the range including equality at the endpoints.
         *
         * This function works on types `T` with equality comparison operators
         * defined.
         *
         * @tparam T The arithmetic type to work with.
         * @param value The value to check for containment in the range.
         * @param min The included start of the range.
         * @param max The included end of the range.
         * @return Whether `value` was within the range.
         */
        template <typename T>
        [[nodiscard]] constexpr bool inRangeClosed(T value, T min, T max)
        {
            return min <= value && value <= max;
        }

        /**
         * @brief Given a value and min/max of an arithemtic type `T`, checks if
         * the value is within the range excluding equality at the endpoints.
         *
         * This function works on types `T` with strict comparison operators
         * defined.
         *
         * @tparam T The arithmetic type to work with.
         * @param value The value to check for containment in the range.
         * @param min the excluded start of the range.
         * @param max the excluded end of the range.
         * @return Whether `value` was within the range.
         */
        template <typename T>
        [[nodiscard]] constexpr bool inRangeOpen(T value, T min, T max)
        {
            return min < value && value < max;
        }
    }
}
