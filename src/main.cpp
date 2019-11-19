
#define mth_ELIMINATION
#define mth_ROW_MAJOR

#include <gtest/gtest.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include <mth/comp.hpp>
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
    ASSERT_EQ(acz[0], 'a') << "vector { 'a', 'c', 'z' }[0] was not 'a'";
    ASSERT_EQ(acz.get(2), 'z') << "vector { 'a', 'c', 'z' }.get(2) was not 'z'";
    ASSERT_EQ(acz.component<0>(), acz.get(0))
        << "mth::tvec::component inconsistent with mth::tvec::get";
}

TEST(Vectors, OutOfRangeThrows)
{
    mth::ivec2 xy;
    ASSERT_THROW(xy.get(4), std::out_of_range);
    ASSERT_THROW((2 * xy + xy).get(7), std::out_of_range);
}

TEST(Vectors, IntVecDefaultInitIsZero)
{
    mth::tvec<int, 7> vector;

    for (size_t i = 0; i < 7; i++)
    {
        ASSERT_EQ(vector[i], 0)
            << "component " << i
            << " of default initialized mth::tvec<int, 7> was not 0";
    }
}

TEST(Vectors, TypeAliasesExist)
{
    static_assert(
        std::is_same_v<mth::dvec2, mth::tvec<double, 2>>,
        "expected alias mth::dvec2 to be mth::tvec<double, 2>");
    static_assert(
        std::is_same_v<mth::ivec4, mth::tvec<int, 4>>,
        "expected alias mth::ivec4 to be mth::tvec<int, 4>");
    static_assert(
        std::is_same_v<mth::fvec3, mth::tvec<float, 3>>,
        "expected alias mth::fvec3 to be mth::tvec<float, 3>");
    static_assert(
        std::is_same_v<mth::uvec2, mth::tvec<unsigned int, 2>>,
        "expected alias mth::uvec2 to be mth::tvec<unsigned int, 2>");

    static_assert(
        std::is_same_v<mth::vec3, mth::fvec3>,
        "expected alias mth::vec3 to be mth::fvec3");
    static_assert(
        std::is_same_v<mth::vec4, mth::fvec4>,
        "expected alias mth::vec4 to be mth::fvec4");
}

TEST(Vectors, AdditionIsComponentWise)
{
    mth::vec3 vectorLHS = {0.1f, 0.5f, -1.f};
    mth::vec3 vectorRHS = {0.2f, -2.3f, 5.9f};

    auto sum = vectorLHS + vectorRHS;

    for (size_t i = 0; i < 3; i++)
    {
        ASSERT_EQ(sum[i], vectorLHS[i] + vectorRHS[i])
            << "float component " << i
            << " of a vector sum was not equal to the sum of the respective "
               "components of the operands";
    }
}

TEST(Vectors, NegationIsComponentWise)
{
    mth::ivec3 vector = {-2, 3, 1};
    auto negation     = -vector;

    for (size_t i = 0; i < 3; i++)
    {
        ASSERT_EQ(negation[i], -vector[i])
            << "component " << i
            << " of negation of a vector was not the negation of it's "
               "corresponding component";
    }
}

TEST(Vectors, SubtractionIsComponentWise)
{
    mth::dvec3 vectorLHS = {-0.9, 0.2, 1.3};
    mth::dvec3 vectorRHS = {0.65, 0.3, -1.9};

    auto difference = vectorLHS - vectorRHS;

    for (size_t i = 0; i < 3; i++)
    {
        ASSERT_EQ(difference[i], vectorLHS[i] - vectorRHS[i])
            << "float component " << i
            << " of a vector difference was not equal to the difference of the "
               "respective "
               "components of the operands";
    }
}

TEST(Vectors, ScalarMultiplicationDistributesComponentWise)
{
    mth::tvec<int, 7> vector = {1, 5, -2, 0, 0, 3, -4};
    int scalar               = 5;

    auto product = scalar * vector;

    static_assert(
        std::is_same_v<decltype(scalar * vector), decltype(vector * scalar)>,
        "expected scalar product to allow commutative syntax");

    for (size_t i = 0; i < 7; i++)
    {
        ASSERT_EQ(product[i], scalar * vector[i])
            << "component " << i
            << " of product of scalar with vector was not the same as the "
               "scalar multiplied by the component";
    }
}

TEST(Vectors, ScalarDivisionDistributesComponentWise)
{
    mth::tvec<float, 4> vector = {0.1f, 0.5f, -0.7f, 1.95f};
    float scalar               = 5.0f;

    auto reduction = vector / scalar;

    for (size_t i = 0; i < 4; i++)
    {
        ASSERT_EQ(reduction[i], vector[i] / scalar)
            << "component " << i
            << " of vector reduced by scalar was not the same as the original "
               "component divided by the scalar";
    }
}

TEST(Vectors, IteratorsCoverComponents)
{
    mth::uvec4 vector = {0U, 5U, 2U, 7U};

    size_t i = 0;
    for (auto elem : vector)
    {
        ASSERT_EQ(elem, vector[i])
            << "iteration " << i
            << " over vector was not equal to the component at the same index";
        i++;
    }
}

TEST(Vectors, HaveComponentAliases)
{
    mth::vec3 v = {0.1f, 4.7f, -1.1f};

    ASSERT_EQ(v.x(), v.get(0))
        << "x value of mth::vec3 was different from component 0";
    ASSERT_EQ(v.y(), v.get(1))
        << "y value of mth::vec3 was different from component 1";
    ASSERT_EQ(v.z(), v.get(2))
        << "z value of mth::vec3 was different from component 2";

    mth::uvec4 rgba = {255U, 128U, 50U, 255U};

    ASSERT_EQ(rgba.r(), 255)
        << "r value of rgba vetor was different from first component";
    ASSERT_EQ(rgba.g(), 128)
        << "g value of rgba vetor was different from second component";
    ASSERT_EQ(rgba.b(), 50)
        << "b value of rgba vetor was different from third component";
    ASSERT_EQ(rgba.a(), 255)
        << "a value of rgba vetor was different from last component";
}

TEST(Vectors, HaveEquality)
{
    mth::tvec<unsigned int, 5> a = {1U, 13U, 298U, 4U, 37U};
    auto b                       = a;
    mth::tvec<unsigned int, 5> c = {3U, 13U, 5U, 49929U, 12U};

    ASSERT_EQ(c, c) << "equality of vectors failed to be reflexive";
    ASSERT_EQ(a, b)
        << "equality of vectors did not respect copy initialization";
    ASSERT_EQ(b, a) << "equality of vectors failed to be symmetric";
    ASSERT_NE(a, c) << "unequal vectors falsely flagged equal";
}

TEST(Vectors, AreMappable)
{
    mth::uvec4 color255 = {183U, 253U, 86U, 255U};
    auto color01        = color255.map([](auto v) {
        return v / 255.0f;
    });

    ASSERT_EQ(color01.get(2), color255.get(2) / 255.0f)
        << "component 2 of the mapped color was not the same as applying the "
           "mapping to component 2 of the pre image";
}

TEST(Vectors, AreMultiMappable)
{
    mth::vec3 coeffs                = {1.0f, 2.0f, -1.5f};
    mth::tvec<mth::vec3, 3> vectors = {mth::vec3{0.1f, 0.2f, 0.3f},
                                       mth::vec3{1.0f, -1.0f, 0.0f},
                                       mth::vec3{0.1f, 0.5f, 0.3f}};
    auto result                     = mth::vec::map(
        [](auto coeff, auto vector) {
            return coeff * vector;
        },
        coeffs,
        vectors);
    ASSERT_EQ(result[0], coeffs[0] * vectors[0])
        << "first component of mth::vec::map result was miscalculated";
}

// TODO: dot product, magnitude

TEST(Complex, CanGetComponents)
{
    mth::comp z{2.0f, 3.0f};

    ASSERT_EQ(z.real(), 2.0f)
        << "extracted real part didn't match value passed to constructor";
    ASSERT_EQ(z.imag(), 3.0f)
        << "extracted real part didn't match value passed to constructor";
}

TEST(Complex, DefaultInitToZero)
{
    mth::icomp z;

    ASSERT_EQ(z.real(), 0)
        << "default int complex number had non-zero real part";
    ASSERT_EQ(z.imag(), 0)
        << "default int complex number had non-zero imaginary part";
}

TEST(Complex, TypeAliasesExist)
{
    static_assert(
        std::is_same_v<mth::dcomp, mth::tcomp<double>>,
        "expected alias mth::dcomp to be mth::tcomp<double>");
    static_assert(
        std::is_same_v<mth::fcomp, mth::tcomp<float>>,
        "expected alias mth::fcomp to be mth::tcomp<float>");
    static_assert(
        std::is_same_v<mth::icomp, mth::tcomp<int>>,
        "expected alias mth::icomp to be mth::tcomp<int>");
    static_assert(
        std::is_same_v<mth::ucomp, mth::tcomp<unsigned int>>,
        "expected alias mth::ucomp to be mth::tcomp<unsigned int>");
    static_assert(
        std::is_same_v<mth::comp, mth::fcomp>,
        "expected alias mth::comp to be mth::fcomp");
}

TEST(Complex, AdditionIsCorrect)
{
    mth::icomp LHS = {1, 3};
    mth::icomp RHS = {-2, 5};

    auto sum = LHS + RHS;

    ASSERT_EQ(sum.real(), LHS.real() + RHS.real())
        << "real part of complex sum was not the sum of the individual real "
           "parts";
    ASSERT_EQ(sum.imag(), LHS.imag() + RHS.imag())
        << "imaginary part of complex sum was not the sum of the individual "
           "imaginary parts";
}

TEST(Complex, NegationIsCorrect)
{
    mth::ucomp z  = {5, 0};
    auto negation = -z;

    ASSERT_EQ(negation.real(), -z.real())
        << "real part of complex negation was not the negation of the real "
           "part";
    ASSERT_EQ(negation.imag(), -z.imag())
        << "imaginary part of complex negation was not the negation of the "
           "imaginary part";
}

TEST(Complex, ConjugationIsCorrect)
{
    mth::dcomp z = {1.24, 2.7899};
    auto conj    = z.conj();

    ASSERT_EQ(conj.real(), z.real()) << "conjugation modified the real part";
    ASSERT_EQ(conj.imag(), -z.imag())
        << "conjugation did not negate the imaginary part";
}

TEST(Complex, SubtractionIsCorrect)
{
    mth::icomp a = {4, -7};
    mth::icomp b = {-1, 0};

    auto diff = a - b;

    ASSERT_EQ(diff.real(), a.real() - b.real())
        << "real part of complex difference was not the difference of the "
           "individual real parts";
    ASSERT_EQ(diff.imag(), a.imag() - b.imag())
        << "imaginary part of complex difference was not the difference of the "
           "individual imaginary parts";
}

TEST(Complex, MultiplicationIsCorrect)
{
    mth::ucomp a = {3, 9};
    mth::ucomp b = {2, 6};

    auto prod = a * b;

    ASSERT_EQ(prod.real(), a.real() * b.real() - a.imag() * b.imag())
        << "real part of complex product was incorrect";
    ASSERT_EQ(prod.imag(), a.real() * b.imag() + a.imag() * b.real())
        << "imaginary part of complex product was incorrect";
}

TEST(Complex, CanFindModulus)
{
    mth::dcomp z = {1.2, 2.9};

    ASSERT_EQ(z.magnSqr(), 1.2 * 1.2 + 2.9 * 2.9)
        << "square modulus of a complex number was incorrect";
    ASSERT_EQ(z.magn(), std::sqrt(z.magnSqr()))
        << "modulus of a complex number was incorrect";
}

TEST(Complex, HaveEquality)
{
    mth::icomp a = {13, 27};
    mth::icomp b = {-9, -1};
    mth::icomp c = a;

    ASSERT_EQ(b, b) << "equality of complex numbers failed to be reflexive";
    ASSERT_EQ(a, c)
        << "equality of complex numbers didn't respect copy initialization";
    ASSERT_NE(a, b) << "unequal complex numbers falsely flagged equal";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
