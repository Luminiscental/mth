
#define m_PRECISION 16

#include <iostream>

#include <m/m.h>

#include <m/vec.h>
#include <m/mat.h>
#include <m/quat.h>
#include <m/polynomial.h>
#include <m/series.h>
#include <m/numeric.h>

size_t factorial(size_t n) {

    if (n == 0) return 1;

    auto result = n;

    for (size_t i = 2; i < n; i++) {

        result *= i;
    }

    return result;
}

auto innerProd3(m::vec3 a, m::vec3 b) {

    m::mat3x1 aMatT(a);
    m::mat1x3 bMat(b);

    m::mat1x1 result = aMatT * bMat;

    return result.get();
}

auto outerProd3(m::vec3 a, m::vec3 b) {

    m::mat1x3 aMat(a);
    m::mat3x1 bMatT(b);

    m::mat3x3 result = aMat * bMatT;

    return result;
}

#define printThing(thing) std::cout << std::endl << thing << std::endl

// TODO: Actual tests

int main() {

    {

        printThing("-- testing complex math functions --");

        auto a = 1.0 + 2.0 * m::i<double>;
        auto b = 2.0 - 3.0 * m::i<double>;

        printThing(a << " * " << b << " = " << (a * b));
        printThing("sqrt(" << a << ") = " << std::sqrt(a));
        printThing("exp(" << b << ") = " << std::exp(b));
        printThing("cos(" << b << ") = " << std::cos(b));

        auto aPolar = a.asPolar();

        printThing("polar form of " << a << " is " << aPolar);
        printThing(aPolar.x() << " * e^(i * " << aPolar.y() << ") = " << (aPolar.x() * std::exp(m::i<double> * aPolar.y())));

        auto c = 1.0 + m::i<double>;
        auto cosC = std::cos(c);
        auto sinC = std::sin(c);

        printThing("c = " << c);
        printThing("cos(c) = " << cosC);
        printThing("sin(c) = " << sinC);

        printThing(cosC << " + i * " << sinC << " = " << (cosC + m::i<double> * sinC));
        printThing("exp(ic) = " << std::exp(m::i<double> * c));
    }

    {

        printThing("-- testing non-square matrix multiplication and vector functions --");

        m::vec3 someVector(1, 3, 5);
        m::vec3 anotherVector = 7.0 * someVector.unit();
        auto inner = innerProd3(someVector, anotherVector);
        auto outer = outerProd3(someVector, anotherVector);

        printThing("inner(" << someVector << ", " << anotherVector << ") = " << inner);
        printThing("outer(" << someVector << ", " << anotherVector << ") =");
        printThing(outer);
    }

    {

        printThing("-- testing matrix / quaternion transformations of 3-vectors --");

        m::vec3 initialPosition(1, 2, 1);

        m::mat4 transformation = m::mat::rotation(m::PI<double> / 3, m::X_AXIS<double>)
                               * m::mat::translation(m::vec3(1.0, 13.0, -2.0));
        m::vec4 positionHandle(initialPosition.x(), initialPosition.y(), initialPosition.z(), 1);

        m::vec3 transformedPosition = (transformation * positionHandle).xyz();

        m::quat additionalRotation = m::quat::rotation(m::TAU<double> / 12, m::vec3(0.5f, 0.5f, -0.1f));
        additionalRotation *= m::quat(2, 3, 1, 5);
        additionalRotation.real() -= 4;

        transformedPosition = additionalRotation.rotate(transformedPosition);

        printThing("Transformation matrix:");
        printThing(m::mat::rotation(additionalRotation) * transformation);

        printThing("Went from " << initialPosition << " to " << transformedPosition);
    }

    {

        printThing("-- testing inverse of large matrix --");

        m::tmat<double, 8, 8> bigMatrix(0.0, 3.0, -1.0, 5.0, 7.0, 6.0, 2.0, -3.0,
                                        1.0, 5.0, 2.0, -3.0, 2.0, -2.0, -2.0, 5.0,
                                        7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,
                                        3.0, 2.0, -1.0, 3.0, -5.0, -6.0, -1.0, 1.0,
                                        5.0, 3.0, 4.0, 3.0, 5.0, -2.0, -3.0, 1.0,
                                        6.0, 4.0, 3.0, 2.0, 7.0, -5.0, -7.0, 11.0,
                                        0.0, 2.0, 0.0, 1.5, 2.5, -3.0, 5.0, 4.0,
                                        2.0, 4.0, 8.0, 12.0, -11.0, 2.0, 5.0, -6.0);

        printThing("A :");
        printThing(bigMatrix);

        printThing("A^{-1} :");
        printThing(bigMatrix.inverse());
    }

    {

        printThing("-- testing matrix edge-cases --");

        m::tmat<double, 2, 3> rectMat(1, 2,
                                     3, 4,
                                     5, 6);

        auto rectProd = rectMat * m::vec2(-1, -2);

        printThing("B :");
        printThing(rectMat);

        printThing("B . (-1, -2) = " << rectProd);

        auto tinyMatrix = m::mat1(3);
        auto tinyInverse = tinyMatrix.inverse();

        printThing("C :");
        printThing(tinyMatrix);

        printThing("C^{-1} :");
        printThing(tinyInverse);

        auto complicatedMatrix = m::cmat4x2(1, 3.0 * m::i<double>, -m::i<double>, 2.0 * m::i<double>,
                                            -2.0 * m::i<double>, 3, 2, 1);

        auto complexVector = m::cvec4(1, m::i<double>, -1, -m::i<double>);

        auto complicatedProd = complicatedMatrix * complexVector;

        printThing("D :");
        printThing(complicatedMatrix);

        printThing("D . " << complexVector << " = " << complicatedProd);
    }

    {

        printThing("-- testing polynomial arithmetic --");

        auto quadratic = m::Polynomial::fromCoeffs(1, 2, 3); // 1 + 2z + 3z^2

        auto roots = quadratic.solve();
        auto degree = quadratic.getDegree();

        printThing("P(z) = " << quadratic);
        printThing("P has degree " << degree);
        printThing("P(z) = 0 -> " << roots);

        auto quintic = m::Polynomial::fromCoeffs(1, 0, 0, 0, 0, 1); // 1 + z^5
        auto prod = quadratic * quintic;

        printThing("Q(z) = " << quintic);
        printThing("P(z) * Q(z) = R(z) = " << prod);

        printThing("R(1) = " << prod.value(1));
        printThing("R has degree " << prod.getDegree());
    }

    {

        printThing("-- testing analytic calculus --");

        auto thing = m::Polynomial::fromCoeffs(1, 2, 1, 2, 1, 2);

        auto d = m::differentiate(thing);
        auto p = m::integrate(thing);

        printThing("A(z) = " << thing);
        printThing("A'(z) = " << d);
        printThing("B'(z) = A(z) and B(0) = 0 -> B(z) = " << p);
    }

    {

        printThing("-- testing limits --");

        auto sincFunc = [] (m::comp z) { return std::sin(z) / z; };
        auto sincLimit = m::numeric::limit(sincFunc, 0);

        printThing("lim z->0 { sin(z)/z } = " << sincLimit);

        auto num = m::Polynomial::fromCoeffs(1, 0, 1);
        auto denom = m::Polynomial::fromCoeffs(0, 1, 2);

        auto rational = [num, denom] (size_t n) {

            return num.value(m::comp(n)) / denom.value(m::comp(n));
        };

        auto rationalFunc = [num, denom] (m::comp z) {

            return num.value(z) / denom.value(z);
        };

        // Providing values between makes it more accurate

        auto lim = m::numeric::limit(rational);
        auto limFunc = m::numeric::limitInfPos(rationalFunc);
        printThing("lim n->infinity (" << num << ") / (" << denom << ") = " << lim);
        printThing("lim z->infinity (" << num << ") / (" << denom << ") = " << limFunc);
    }

    {

        printThing("-- testing polynomial interpolation --");

        std::vector<m::cvec2> points;

        points.push_back(m::cvec2(0, 1));
        points.push_back(m::cvec2(1, 2));
        points.push_back(m::cvec2(2, 1));

        auto pol = m::Polynomial::interpolate(points);

        printThing("(0, 1) and (1, 2) and (2, 1) gives " << pol);

        std::vector<m::cvec2> line;

        line.push_back(m::cvec2(0.5, 0.1));
        line.push_back(m::cvec2(0.25, 0.1));
        line.push_back(m::cvec2(0.125, 0.1));
        line.push_back(m::cvec2(0.05, 0.1));

        auto shouldBeConst = m::Polynomial::interpolate(line);

        printThing("Points where y=0.1 gives " << shouldBeConst);
    }

    {

        printThing("-- testing series --");

        auto expGen = [] (size_t n) {

            return m::comp(1) / m::comp(factorial(n));
        };

        m::Series expSeries(expGen);
        printThing("e is roughly " << expSeries.getValue(m::comp(1)));

        auto arctanCoeffs = [] (size_t n) {

            if (n % 2 == 0) return m::comp(0);

            m::comp sign = 1;

            if (n % 4 == 3) sign = -1;

            return sign / m::comp(n);
        };

        auto machinCoeffs = [arctanCoeffs] (size_t n) {

            if (n % 2 == 0) return m::comp(0);

            auto orig = arctanCoeffs(n);
            auto over5 = orig / std::pow(m::comp(5), m::comp(n));
            auto over239 = orig / std::pow(m::comp(239), m::comp(n));

            return m::comp(4) * (m::comp(4) * over5 - over239);
        };

        m::Series machin(machinCoeffs);
        printThing("pi is roughly " << machin.getValue(m::comp(1)));
    }

    std::cout << std::endl;

    return 0;
}

#undef printThing
