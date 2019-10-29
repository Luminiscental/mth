
#define mth_ELIMINATION
#define mth_ROW_MAJOR

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>

#include <mth/mth.hpp>
#include <mth/vec.hpp>

TEST(GenericUtilities, CanGetIntZero)
{
    auto zero = mth::util::zero<int>();

    static_assert(
        std::is_same_v<decltype(zero), int>,
        "expected mth::util::zero<int>() to be a int");
    ASSERT_EQ(zero, 0) << "mth::util::zero<int>() is not 0";
}

TEST(GenericUtilities, CanGetDoubleZero)
{
    auto zero = mth::util::zero<double>();

    static_assert(
        std::is_same_v<decltype(zero), double>,
        "expected mth::util::zero<double>() to be a double");
    ASSERT_EQ(zero, 0.0) << "mth::util::zero<double>() is not 0.0";
}

TEST(GenericUtilities, IntZeroIsZero)
{
    auto zero = mth::util::zero<int>();
    ASSERT_TRUE(mth::util::isZero(zero))
        << "mth::util::zero<int>() doesn't satisfy mth::util::isZero()";
}

TEST(GenericUtilities, DoubleZeroIsZero)
{
    auto zero = mth::util::zero<double>();
    ASSERT_TRUE(mth::util::isZero(zero))
        << "mth::util::zero<double>() doesn't satisfy mth::util::isZero()";
}

TEST(GenericUtilities, IntOneIsInZeroToTwoOpen)
{
    ASSERT_TRUE(mth::util::inRangeOpen<int>(1, 0, 2))
        << "mth::util::inRangeOpen<int> thinks 1 isn't in (0,2)";
}

TEST(GenericUtilities, IntOneIsInOneToThreeClosed)
{
    ASSERT_TRUE(mth::util::inRangeClosed<int>(1, 1, 3))
        << "mth::util::inRangeClosed<int> thinks 1 isn't in [1,3]";
}

TEST(GenericUtilities, DoubleTenthInZeroToOneOpen)
{
    ASSERT_TRUE(mth::util::inRangeOpen<double>(0.1, 0, 1))
        << "mth::util::inRangeOpen<double> thinks 0.1 isn't in (0,1)";
}

TEST(Vectors, CanGetComponent)
{
    mth::tvec<char, 3> acz = {'a', 'c', 'z'};

    ASSERT_EQ(acz.get(1), 'c') << "vector { 'a', 'c', 'z' }.get(1) was not 'c'";
}

TEST(Vectors, IntVecDefaultInitIsZero)
{
    mth::tvec<int, 7> vector;

    for (size_t i = 0; i < 7; i++)
    {
        ASSERT_EQ(vector.get(i), 0)
            << "component " << i
            << " of default initialized mth::tvec<int, 7> was not 0";
    }
}

TEST(Vectors, AdditionIsComponentWise)
{
    mth::tvec<float, 3> vectorLHS = {0.1f, 0.5f, -1.f};
    mth::tvec<float, 3> vectorRHS = {0.2f, -2.3f, 5.9f};

    auto sum = vectorLHS + vectorRHS;

    static_assert(
        std::is_same_v<decltype(sum), mth::tvec<float, 3>>,
        "expected sum of two vectors to be of the same type as the operands");

    for (size_t i = 0; i < 3; i++)
    {
        ASSERT_EQ(sum.get(i), vectorLHS.get(i) + vectorRHS.get(i))
            << "float component " << i
            << " of a vector sum was not equal to the sum of the respective "
               "components of the operands";
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
