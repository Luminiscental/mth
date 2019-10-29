#pragma once

namespace mth
{
    namespace util
    {
        // Given an arithmetic type T, tries to give the representation
        // of zero in that type.
        template <typename T>
        [[nodiscard]] constexpr T zero() noexcept
        {
            return T{} * 0;
        }
    }
}
