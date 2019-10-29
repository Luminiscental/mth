
#define mth_ELIMINATION
#define mth_ROW_MAJOR

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>

#include <mth/mth.hpp>

TEST(GenericUtilities, CanGetDoubleZero)
{
    auto zero = mth::util::zero<double>();

    static_assert(
        std::is_same_v<decltype(zero), double>,
        "expected mth::util::zero<double>() to be a double");
    ASSERT_EQ(zero, 0.0) << "mth::util::zero<double>() is not 0.0";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
