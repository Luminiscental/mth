#ifndef __m_polynomial_h__
#define __m_polynomial_h__

#include <unordered_set>
#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <m/comp.h>

// TODO: Move to comp and add for other classes
namespace std {

    template<>
    struct hash<m::comp> {

        size_t operator()(const m::comp &x) const;
    };
}

namespace m {

    class ComplexSolutions;

    class PolynomialDegree;

    class Polynomial;

    auto &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

    auto operator==(const PolynomialDegree &lhs, const PolynomialDegree &rhs);
    auto operator!=(const PolynomialDegree &lhs, const PolynomialDegree &rhs);

    auto &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

    auto operator+(const Polynomial &lhs, const comp &rhs);
    auto operator+(const comp &lhs, const Polynomial &rhs);
    auto operator+(const Polynomial &lhs, const Polynomial &rhs);

    auto operator-(const Polynomial &rhs);
    auto operator-(const Polynomial &lhs, const comp &rhs);
    auto operator-(const comp &lhs, const Polynomial &rhs);
    auto operator-(const Polynomial &lhs, const Polynomial &rhs);

    auto operator*(const Polynomial &lhs, const comp &rhs);
    auto operator*(const comp &lhs, const Polynomial &rhs);
    auto operator*(const Polynomial &lhs, const Polynomial &rhs);

    auto operator/(const Polynomial &lhs, const comp &rhs);

    auto operator==(const Polynomial &lhs, const Polynomial &rhs);
    auto operator!=(const Polynomial &lhs, const Polynomial &rhs);

    auto &operator<<(std::ostream &lhs, const Polynomial &rhs);

    class ComplexSolutions {

    private:

        inline ComplexSolutions() noexcept;
        inline ComplexSolutions(std::unordered_set<comp> finiteSet) noexcept;

    public:

        std::unordered_set<comp> solutionSet;
        bool inf;

        inline static ComplexSolutions empty() noexcept;

        inline static ComplexSolutions finite(std::unordered_set<comp> finiteSet) noexcept;
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept;

        inline static ComplexSolutions infinite() noexcept;

        friend auto &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

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

        friend auto &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

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

        auto getCoeffs();
        auto getCoeffs() const;
        auto getDegree();
        auto getDegree() const;

        comp value(comp x) const;
        auto solve();

        auto getCoeff(size_t index) const;
        void setCoeff(size_t index, const comp &value);

        // TODO: Unfriend and use getters/setter (then also in other classes)

        friend auto operator+(const Polynomial &lhs, const comp &rhs);
        friend auto operator+(const comp &lhs, const Polynomial &rhs);
        friend auto operator+(const Polynomial &lhs, const Polynomial &rhs);

        friend auto operator-(const Polynomial &rhs);
        friend auto operator-(const Polynomial &lhs, const comp &rhs);
        friend auto operator-(const comp &lhs, const Polynomial &rhs);
        friend auto operator-(const Polynomial &lhs, const Polynomial &rhs);

        friend auto operator*(const Polynomial &lhs, const comp &rhs);
        friend auto operator*(const comp &lhs, const Polynomial &rhs);
        friend auto operator*(const Polynomial &lhs, const Polynomial &rhs);

        friend auto operator/(const Polynomial &lhs, const comp &rhs);

        friend auto operator==(const Polynomial &lhs, const Polynomial &rhs);
        friend auto operator!=(const Polynomial &lhs, const Polynomial &rhs);

        friend auto &operator<<(std::ostream &lhs, const Polynomial &rhs);
    };

    Polynomial differentiate(const Polynomial &polynomial);
    Polynomial integrate(const Polynomial &polynomial);
}

#include <m/polynomial_impl.h>

#endif
