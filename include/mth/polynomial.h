#ifndef mth_polynomial_h__
#define mth_polynomial_h__

/* <mth/polynomial.h> - polynomial header
 *      This includes the class Polynomial which stores a polynomial with
 *      complex coefficients. Ring arithmetic operators are defined as well
 *      as member functions to find values such as the degree and roots of the
 *      stored polynomial. The polynomial can be evaluated at a point with
 *      value() or converted to a complex function represented with
 *      std::function. Differentiation and integration is also implemented.
 * 
 *      This header also includes classes to represent a solution set of
 *      complex numbers and the degree of a polynomial, including the
 *      infinite / empty degenerate cases.
 */

#include <unordered_set>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <mth/mth.h>
#include <mth/vec.h>
#include <mth/comp.h>

namespace mth {

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

    // Class representing a set of complex numbers that solve an equation
    class ComplexSolutions {

    private:

        char variableName = 'z';

        std::unordered_set<comp> solutionSet;
        bool inf = false;

        // Initializers private to hide implementation
        ComplexSolutions() = default;
        ComplexSolutions(std::unordered_set<comp> finiteSet) noexcept;

    public:

        // Check whether a number is within this set
        bool contains(const comp &z) const;

        // Check whether this set is infinite
        bool isInfinite() const;

        // Create as empty
        static ComplexSolutions empty() noexcept;

        // Create from a finite set of values as an unordered_set
        static ComplexSolutions finite(std::unordered_set<comp> finiteSet) noexcept;
        
        // Create from a finite set of values as variadic arguments
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept {

            return ComplexSolutions(std::unordered_set<mth::comp>({args...}));
        }

        // Create as infinite
        static ComplexSolutions infinite() noexcept;

        // Overwrite the variable name
        ComplexSolutions setVariableName(char newName);

        // Get the variable name
        char getVariableName() const;

        friend std::ostream &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

        friend class Polynomial;
    };

    // Class to represent the degree of a polynomial
    class PolynomialDegree {

    private:

        size_t value;
        bool inf = false;

    public:

        // Initialize as a finite value
        PolynomialDegree(size_t value);

        // Check if the degree is infinite
        bool isInfinite() const;

        // Get as a finite value (0 when infinite)
        size_t getValue() const;

        // Static factory for infinite degree to hide implementation
        static PolynomialDegree infinite();

        // Friend operators to allow access to inf
        
        friend bool operator==(const PolynomialDegree &lhs, const PolynomialDegree &rhs);
        friend bool operator!=(const PolynomialDegree &lhs, const PolynomialDegree &rhs);

        friend std::ostream &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

        friend class Polynomial;
    };

    class Polynomial {

    private:

        // TODO: Use a string with more sophisticated validation
        char variableName = 'z';

        std::vector<comp> coeffs;

        ComplexSolutions roots = ComplexSolutions::empty();
        bool rootsValid = true;

        PolynomialDegree degree = PolynomialDegree::infinite();
        bool degreeValid = true;

        // Trims zero coefficients from the back of coeffs
        void updateValues();

        // Calculate the degree by iterating to the last non-zero coefficient
        void updateDegree() noexcept;

    public:

        // Initialize to 0
        Polynomial() noexcept {}

        // Create from a list of coefficients
        static Polynomial fromCoeffs(std::vector<comp> coeffs) noexcept;

        // Create from coefficients as a variadic argument list
        template <typename ...Q>
        static Polynomial fromCoeffs(Q... args) {

            return fromCoeffs(std::vector<comp> {static_cast<comp>(args)...});
        }

        // Convert to a complex function
        operator std::function<comp(comp)>() const;

        // Call as a complex function
        comp operator()(const comp &z) const;

        // Getters for the vector of coefficients
        const std::vector<comp> &getCoeffs();
        const std::vector<comp> &getCoeffs() const;

        // Recalculates the degree if needed
        PolynomialDegree getDegree();
        PolynomialDegree getDegree() const;

        // Evaluate at a point
        comp value(comp z) const;

        // Solve setting the polynomial to zero
        ComplexSolutions solve();

        // Get the coefficient at an index
        comp getCoeff(size_t index) const;

        // Overwrite the coefficient at an index
        void setCoeff(size_t index, const comp &value);

        // Overwrite the variable name
        Polynomial setVariableName(char newName);

        // Get the variable name
        char getVariableName() const;

        // Returns the lagrange interpolation polynomial intersecting the specified points in range
        static Polynomial interpolate(const std::vector<cvec2> &points, size_t first, size_t last);

        // Overload to default first=0, last=size()-1
        static Polynomial interpolate(const std::vector<cvec2> &points);

        // Compound assignment operators

        Polynomial &operator+=(const Polynomial &rhs);
        Polynomial &operator-=(const Polynomial &rhs);
        Polynomial &operator*=(const Polynomial &rhs);

        Polynomial &operator+=(const comp &rhs);
        Polynomial &operator-=(const comp &rhs);
        Polynomial &operator*=(const comp &rhs);
        Polynomial &operator/=(const comp &rhs);

        // Friend operators

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

    // Differentiation and integration of polynomials

    Polynomial differentiate(const Polynomial &polynomial);
    Polynomial integrate(const Polynomial &polynomial);
}

#endif
