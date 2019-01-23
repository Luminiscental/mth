#ifndef __m_polynomial_h__
#define __m_polynomial_h__

/* <m/polynomial.h> - polynomial header
 *      This include the class Polynomial which stores a polynomial with complex coefficients. Addition, subtraction
 *      and multiplication operators are defined as well as member functions to find values such as the degree and roots
 *      of the stored polynomial. The polynomial can be evaluated at a point with value() or converted to a complex function
 *      represented with std::function. Differentiation and integration is also implemented.
 * 
 *      This header also includes classes to represent a solution set of complex numbers and the degree of a polynomial,
 *      including the infinite / empty degenerate cases.
 */

#include <unordered_set>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <m/m.h>
#include <m/comp.h>

namespace m {

    // Forward declaration for friending

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

        // Initializers private to hide implementation
        ComplexSolutions() noexcept;
        ComplexSolutions(std::unordered_set<comp> finiteSet) noexcept;

    public:

        // TODO: Use getters instead
        std::unordered_set<comp> solutionSet;
        bool inf;

        // Static factory functions for producting solution sets
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

        // Static factory for infinite degree to hide implementation
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

        // Trims zero coefficients from the back of coeffs
        void updateValues();

        void updateDegree() noexcept;

    public:

        Polynomial() noexcept;

        static Polynomial fromCoeffs(std::vector<comp> coeffs) noexcept;

        template <typename ...Q>
        static Polynomial fromCoeffs(Q... args);

        operator std::function<comp(comp)>() const;
        comp operator()(const comp &z) const;

        std::vector<comp> getCoeffs();
        std::vector<comp> getCoeffs() const;

        // Recalculates the degree if needed
        PolynomialDegree getDegree();
        PolynomialDegree getDegree() const;

        comp value(comp x) const;
        ComplexSolutions solve();

        comp getCoeff(size_t index) const;
        void setCoeff(size_t index, const comp &value);

        // Returns the lagrange interpolation polynomial intersecting the specified points
        static Polynomial interpolate(const std::vector<cvec2> &points);
        static Polynomial interpolate(const std::vector<cvec2> &points, size_t first, size_t last);

        Polynomial &operator+=(const Polynomial &rhs);
        Polynomial &operator-=(const Polynomial &rhs);
        Polynomial &operator*=(const Polynomial &rhs);

        Polynomial &operator+=(const comp &rhs);
        Polynomial &operator-=(const comp &rhs);
        Polynomial &operator*=(const comp &rhs);
        Polynomial &operator/=(const comp &rhs);

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
