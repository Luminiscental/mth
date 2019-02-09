
#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>

#include <m/comp.h>
#include <m/quat.h>

#include <m/vec.h>
#include <m/mat.h>

#include <m/polynomial.h>
#include <m/powerseries.h>
#include <m/numeric.h>

#define m_ASSERT_ZERO(a) ASSERT_TRUE(m::util::isZero(a)) \
    << "Expected " << #a << " which is " << a << " to be zero" << std::endl;

#define m_ASSERT_EQ(a, b) ASSERT_TRUE(m::util::isEqual(a, b)) \
    << "Expected but did not find equality of" << std::endl \
    << "    " << #a << " which is " << a << std::endl \
    << "and" << std::endl \
    << "    " << #b << " which is " << b;

TEST(CompTest, DefaultInitsToZero) {

    m::comp z;

    m_ASSERT_EQ(z.real(), 0.0);
    m_ASSERT_EQ(z.imag(), 0.0);
}

TEST(CompTest, FillsValuesCorrectly) {

    auto z = m::comp::fromCartesian(-1, 2);

    m_ASSERT_EQ(z.real(), -1.0);
    m_ASSERT_EQ(z.imag(), 2.0);
}

TEST(CompTest, ConvertsFromPolarCorrectly) {

    auto diag = m::comp::fromPolar(std::sqrt(2.0), m::PI<double> / 4.0);

    m_ASSERT_EQ(diag.real(), 1.0);
    m_ASSERT_EQ(diag.imag(), 1.0);
}

TEST(CompTest, SumsElementWise) {

    auto p = m::comp::fromCartesian(-7.0, 3.5);
    auto q = m::comp::fromCartesian(3.0, 2.4);

    auto sumPQ = p + q;

    m_ASSERT_EQ(sumPQ, m::comp::fromCartesian(3.0 - 7.0, 3.5 + 2.4));
}

TEST(CompTest, MultipliesCorrectly) {

    auto a = m::comp::fromCartesian(3.0, 2.0);
    auto b = m::comp::fromCartesian(-1.0, 1.0);

    auto prodAB = a * b;

    m_ASSERT_EQ(prodAB, m::comp::fromCartesian(-5.0, 1.0));
}

TEST(CompTest, ConvertsToPolarCorrectly) {

    auto a = m::comp::fromPolar(2.0, 3.0);
    auto b = m::comp::fromPolar(-1.0, 1.1);

    auto prodAB = a * b;

    m_ASSERT_EQ(prodAB, m::comp::fromPolar(-2.0, 4.1));
}

TEST(VecTest, DefaultInitsToZero) {

    m::vec2 new_vec;

    double x = new_vec.x();
    double y = new_vec.y();

    m_ASSERT_EQ(x, 0.0);
    m_ASSERT_EQ(y, 0.0);
}

TEST(VecTest, FillsValuesCorrectly) {

    m::ivec7 seq(1, 2, 3, 4, 5, 6, 7);

    for (int i = 0; i < 7; i++) {

        double elem = seq.get(i);
        double oneIndexed = i + 1;
        m_ASSERT_EQ(elem, oneIndexed);
    }
}

TEST(VecTest, IteratesFully) {

    m::vec5 list_vec(1.2, 1.3, 1.4, 1.5, 1.6);

    double result = 0.0;

    for (auto elem : list_vec) {

        result += elem;
    }

    m_ASSERT_EQ(result, 1.2 + 1.3 + 1.4 + 1.5 + 1.6);
}

TEST(VecTest, DotProdCorrectly) {

    auto x = m::ivec3(1, 3, 2);
    auto y = m::ivec3(2, 4, -1);

    auto dotXY = x.dot(y);

    m_ASSERT_EQ(dotXY, 12);
}

TEST(VecTest, CrossProdCorrectly) {

    m::vec3 p = m::vec3(1.0, 2.0, -1.0);
    m::vec3 q = m::vec3(2.0, 1.0, -2.0);

    m::vec3 crossPQ = m::vec::cross(p, q);

    m_ASSERT_ZERO(crossPQ.dot(p));
    m_ASSERT_ZERO(crossPQ.dot(q));
}

int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

