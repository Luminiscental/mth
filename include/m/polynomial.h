#ifndef __m_polynomial_h__
#define __m_polynomial_h__

#include <unordered_set>
#include <array>
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

    template <size_t N>
    class Polynomial;

    auto &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

    auto &operator<<(std::ostream &lhs, const PolynomialDegree &rhs);

    template <size_t N>
    auto operator+(const Polynomial<N> &lhs, const comp &rhs);

    template <size_t N>
    auto operator+(const comp &lhs, const Polynomial<N> &rhs);

    template <size_t N, size_t M>
    auto operator+(const Polynomial<N> &lhs, const Polynomial<M> &rhs);

    template <size_t N>
    auto operator-(const Polynomial<N> &rhs);

    template <size_t N>
    auto operator-(const Polynomial<N> &lhs, const comp &rhs);

    template <size_t N>
    auto operator-(const comp &lhs, const Polynomial<N> &rhs);

    template <size_t N, size_t M>
    auto operator-(const Polynomial<N> &lhs, const Polynomial<M> &rhs);

    template <size_t N>
    auto operator*(const Polynomial<N> &lhs, const comp &rhs);

    template <size_t N>
    auto operator*(const comp &lhs, const Polynomial<N> &rhs);

    template <size_t N, size_t M>
    auto operator*(const Polynomial<N> &lhs, const Polynomial<M> &rhs);

    template <size_t N>
    auto operator/(const Polynomial<N> &lhs, const comp &rhs);

    template <size_t N, size_t M>
    auto operator==(const Polynomial<N> &lhs, const Polynomial<M> &rhs);

    template <size_t N, size_t M>
    auto operator!=(const Polynomial<N> &lhs, const Polynomial<M> &rhs);

    template <size_t N>
    auto &operator<<(std::ostream &lhs, const Polynomial<N> &rhs);

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

        template <size_t N>
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
    };

    template <size_t N>
    class Polynomial;

    template<>
    class Polynomial<0> {

    private:

        comp coeff;
        ComplexSolutions roots;
        bool rootsValid;

        Polynomial() noexcept;

    public:

        template <size_t M, typename std::enable_if<(M > 0), int>::type = 0>
        Polynomial(std::array<comp, M> coeffs) noexcept;

        Polynomial(comp value) noexcept;

        comp getCoeff() const;

        auto actualDegree() const;

        auto value(comp x);
        auto solve();

        template <size_t M>
        friend auto operator+(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator+(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator+(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator-(const Polynomial<M> &rhs);

        template <size_t M>
        friend auto operator-(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator-(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator-(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator*(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator*(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator*(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator/(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M, size_t O>
        friend auto operator==(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M, size_t O>
        friend auto operator!=(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto &operator<<(std::ostream &lhs, const Polynomial<M> &rhs);
    };

    template <size_t N>
    class Polynomial {

    private:

        std::array<comp, N + 1> coeffs;
        ComplexSolutions roots;
        bool rootsValid;

        Polynomial() noexcept;

    public:


        template <size_t M, typename std::enable_if<(M > N), int>::type = 0>
        Polynomial(std::array<comp, M> coeffs) noexcept;

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type = 0>
        Polynomial(Q... args) noexcept;

        std::array<comp, N + 1> getCoeffs() const;

        auto actualDegree() const;

        auto value(comp x);
        auto solve();

        template <size_t M>
        friend auto operator+(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator+(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator+(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator-(const Polynomial<M> &rhs);

        template <size_t M>
        friend auto operator-(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator-(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator-(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator*(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto operator*(const comp &lhs, const Polynomial<M> &rhs);

        template <size_t M, size_t O>
        friend auto operator*(const Polynomial<M> &lhs, const Polynomial<O> &rhs);

        template <size_t M>
        friend auto operator/(const Polynomial<M> &lhs, const comp &rhs);

        template <size_t M>
        friend auto &operator<<(std::ostream &lhs, const Polynomial<M> &rhs);
    };
}

#include <m/polynomial_impl.h>

#endif
