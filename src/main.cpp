
#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>

#include <mth/comp.h>
#include <mth/quat.h>

#include <mth/vec.h>
#include <mth/mat.h>

#include <mth/polynomial.h>
#include <mth/powerseries.h>
#include <mth/numeric.h>

#define mth_ASSERT_ZERO(a) ASSERT_TRUE(mth::util::isZero(a)) \
    << "Expected " << #a << " which is " << a << " to be zero" << std::endl;

#define mth_ASSERT_EQ(a, b) ASSERT_TRUE(mth::util::isEqual(a, b)) \
    << "Expected but did not find equality of" << std::endl \
    << "    " << #a << " which is " << a << std::endl \
    << "and" << std::endl \
    << "    " << #b << " which is " << b;

TEST(CompTest, DefaultInitsToZero) {

    mth::comp z;

    mth_ASSERT_EQ(z.real(), 0.0);
    mth_ASSERT_EQ(z.imag(), 0.0);
}

TEST(CompTest, FillsValuesCorrectly) {

    auto z = mth::comp::fromCartesian(-1, 2);

    mth_ASSERT_EQ(z.real(), -1.0);
    mth_ASSERT_EQ(z.imag(), 2.0);
}

TEST(CompTest, ConvertsFromPolarCorrectly) {

    auto diag = mth::comp::fromPolar(std::sqrt(2.0), mth::PI<double> / 4.0);

    mth_ASSERT_EQ(diag.real(), 1.0);
    mth_ASSERT_EQ(diag.imag(), 1.0);
}

TEST(CompTest, SumsElementWise) {

    auto p = mth::comp::fromCartesian(-7.0, 3.5);
    auto q = mth::comp::fromCartesian(3.0, 2.4);

    auto sumPQ = p + q;

    mth_ASSERT_EQ(sumPQ, mth::comp::fromCartesian(3.0 - 7.0, 3.5 + 2.4));
}

TEST(CompTest, MultipliesCorrectly) {

    auto a = mth::comp::fromCartesian(3.0, 2.0);
    auto b = mth::comp::fromCartesian(-1.0, 1.0);

    auto prodAB = a * b;

    mth_ASSERT_EQ(prodAB, mth::comp::fromCartesian(-5.0, 1.0));
}

TEST(CompTest, ConvertsToPolarCorrectly) {

    auto a = mth::comp::fromPolar(2.0, 3.0);
    auto b = mth::comp::fromPolar(-1.0, 1.1);

    auto prodAB = a * b;

    mth_ASSERT_EQ(prodAB, mth::comp::fromPolar(-2.0, 4.1));
}

TEST(CompTest, AbsOfZeroIsZero) {

    mth::comp z;

    mth_ASSERT_ZERO(z.abs());
}

TEST(CompTest, ArgConsistentWithPolar) {

    mth::comp w = mth::comp::fromCartesian(1.0, 2.0);

    mth::vec2 asPolar = w.asPolar();
    double arg = w.arg();

    mth_ASSERT_EQ(arg, asPolar.y());
}

TEST(CompTest, AbsConsistentWithPolar) {

    mth::comp w = mth::comp::fromCartesian(2.0, -3.0);

    mth::vec2 asPolar = w.asPolar();
    double abs = w.abs();

    mth_ASSERT_EQ(abs, asPolar.x());
}

TEST(CompTest, InverseIsInverse) {

    mth::comp z = mth::comp::fromCartesian(13.0, -2.0);
    mth::comp zInv = z.inverse();

    mth_ASSERT_EQ(z * zInv, (mth::comp) 1.0);
}

TEST(VecTest, DefaultInitsToZero) {

    mth::vec2 new_vec;

    double x = new_vec.x();
    double y = new_vec.y();

    mth_ASSERT_EQ(x, 0.0);
    mth_ASSERT_EQ(y, 0.0);
}

TEST(VecTest, FillsValuesCorrectly) {

    mth::ivec7 seq(1, 2, 3, 4, 5, 6, 7);

    for (int i = 0; i < 7; i++) {

        double elem = seq.get(i);
        double oneIndexed = i + 1;
        mth_ASSERT_EQ(elem, oneIndexed);
    }
}

TEST(VecTest, IteratesFully) {

    mth::vec5 list_vec(1.2, 1.3, 1.4, 1.5, 1.6);

    double result = 0.0;

    for (auto elem : list_vec) {

        result += elem;
    }

    mth_ASSERT_EQ(result, 1.2 + 1.3 + 1.4 + 1.5 + 1.6);
}

TEST(VecTest, DotProdCorrectly) {

    auto x = mth::ivec3(1, 3, 2);
    auto y = mth::ivec3(2, 4, -1);

    auto dotXY = x.dot(y);

    mth_ASSERT_EQ(dotXY, 12);
}

TEST(VecTest, CrossProdCorrectly) {

    mth::vec3 p = mth::vec3(1.0, 2.0, -1.0);
    mth::vec3 q = mth::vec3(2.0, 1.0, -2.0);

    mth::vec3 crossPQ = mth::vec::cross(p, q);

    mth_ASSERT_ZERO(crossPQ.dot(p));
    mth_ASSERT_ZERO(crossPQ.dot(q));
}

TEST(VecTest, ScalesComponenetWise) {

    mth::vec5 a = mth::vec5(1.0, 2.0, 3.0, 4.0, 5.0);

    mth::vec5 doubled = 2.0 * a;

    mth_ASSERT_EQ(doubled, mth::vec5(2.0, 4.0, 6.0, 8.0, 10.0));
}

TEST(VecTest, SqrMagnMatchesPythag) {

    mth::vec6 b(2.0, -1.0, 13.0, 14.0, -2.5, 1.11);

    double pythag = 0.0;

    for (auto elem : b) {

        pythag += elem * elem;
    }

    mth_ASSERT_EQ(pythag, b.magnSqr());
}

// TODO: Test mat
// TODO: Test quat
// TODO: Test polynomial
// TODO: Test powerseries
// TODO: Test numeric functions

int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

