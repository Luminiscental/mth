#ifndef __m_polynomial_h__
#define __m_polynomial_h__

#include <unordered_set>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <m/comp.h>

namespace m {

    class ComplexSolutions;

    class PolynomialDegree;

    class Polynomial;

    std::ostream &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

    bool operator==(const PolynomialDegree &lhs, const PolynomialDegree &rhs);
    bool operator!=(const PolynomialDegree &lhs, const PolynomialDegree &rhs);

    std::ostream &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

    Polynomial operator+(const Polynomial &lhs, const comp &rhs);
    Polynomial operator+(const comp &lhs, const Polynomial &rhs);
    Polynomial operator+(const Polynomial &lhs, const Polynomial &rhs);

    Polynomial operator-(const Polynomial &rhs);
    Polynomial operator-(const Polynomial &lhs, const comp &rhs);
    Polynomial operator-(const comp &lhs, const Polynomial &rhs);
    Polynomial operator-(const Polynomial &lhs, const Polynomial &rhs);

    Polynomial operator*(const Polynomial &lhs, const comp &rhs);
    Polynomial operator*(const comp &lhs, const Polynomial &rhs);
    Polynomial operator*(const Polynomial &lhs, const Polynomial &rhs);

    Polynomial operator/(const Polynomial &lhs, const comp &rhs);

    bool operator==(const Polynomial &lhs, const Polynomial &rhs);
    bool operator!=(const Polynomial &lhs, const Polynomial &rhs);

    std::ostream &operator<<(std::ostream &lhs, const Polynomial &rhs);

    class ComplexSolutions {

    private:

        ComplexSolutions() noexcept;
        ComplexSolutions(std::unordered_set<comp> finiteSet) noexcept;

    public:

        std::unordered_set<comp> solutionSet;
        bool inf;

        static ComplexSolutions empty() noexcept;

        static ComplexSolutions finite(std::unordered_set<comp> finiteSet) noexcept;
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept;

        static ComplexSolutions infinite() noexcept;

        friend std::ostream &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

        friend class Polynomial;
    };

    class PolynomialDegree {

    private:

        size_t value;
        bool inf;

    public:

        PolynomialDegree(size_t value);

        bool isInfinite() const;
        size_t getValue() const;

        static PolynomialDegree infinite();

        friend bool operator==(const PolynomialDegree &lhs, const PolynomialDegree &rhs);
        friend bool operator!=(const PolynomialDegree &lhs, const PolynomialDegree &rhs);

        friend std::ostream &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

        friend class Polynomial;
    };

    class Polynomial {

    private:

        std::vector<comp> coeffs;

        ComplexSolutions roots;
        bool rootsValid;

        PolynomialDegree degree;
        bool degreeValid;

        void updateValues();
        void updateDegree() noexcept;

    public:

        Polynomial() noexcept;

        Polynomial(std::vector<comp> coeffs) noexcept;

        template <typename ...Q>
        Polynomial(Q... args) noexcept;

        operator std::function<comp(comp)>() const;

        std::vector<comp> getCoeffs();
        std::vector<comp> getCoeffs() const;
        PolynomialDegree getDegree();
        PolynomialDegree getDegree() const;

        comp value(comp x) const;
        ComplexSolutions solve();

        comp getCoeff(size_t index) const;
        void setCoeff(size_t index, const comp &value);

        friend Polynomial operator+(const Polynomial &lhs, const comp &rhs);
        friend Polynomial operator+(const comp &lhs, const Polynomial &rhs);
        friend Polynomial operator+(const Polynomial &lhs, const Polynomial &rhs);

        friend Polynomial operator-(const Polynomial &rhs);
        friend Polynomial operator-(const Polynomial &lhs, const comp &rhs);
        friend Polynomial operator-(const comp &lhs, const Polynomial &rhs);
        friend Polynomial operator-(const Polynomial &lhs, const Polynomial &rhs);

        friend Polynomial operator*(const Polynomial &lhs, const comp &rhs);
        friend Polynomial operator*(const comp &lhs, const Polynomial &rhs);
        friend Polynomial operator*(const Polynomial &lhs, const Polynomial &rhs);

        friend Polynomial operator/(const Polynomial &lhs, const comp &rhs);

        friend bool operator==(const Polynomial &lhs, const Polynomial &rhs);
        friend bool operator!=(const Polynomial &lhs, const Polynomial &rhs);

        friend std::ostream &operator<<(std::ostream &lhs, const Polynomial &rhs);
    };

    Polynomial differentiate(const Polynomial &polynomial);
    Polynomial integrate(const Polynomial &polynomial);
}

#include <m/polynomial_impl.h>

#endif
