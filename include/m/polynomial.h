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

    template <size_t N>
    class Polynomial;

    auto &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

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
    auto &operator<<(std::ostream &lhs, const Polynomial<N> &rhs);

    class ComplexSolutions {

    private:

        std::unordered_set<comp> solutionSet;
        bool inf;

        inline ComplexSolutions() noexcept;
        inline ComplexSolutions(std::unordered_set<comp> finiteSet) noexcept;

    public:

        friend auto &operator<<(std::ostream &lhs, const ComplexSolutions &rhs);

        inline static ComplexSolutions empty() noexcept;

        inline static ComplexSolutions finite(std::unordered_set<comp> finiteSet) noexcept;
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept;

        inline static ComplexSolutions infinite() noexcept;
    };

    template <size_t N>
    class Polynomial;

    template<>
    class Polynomial<0> {

    private:

        comp coeff;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > 0), int>::type = 0>
        Polynomial(std::array<comp, M> coeffs) noexcept;

        Polynomial(comp value) noexcept;

        auto value(comp x);

        auto solve();
    };

    template <size_t N>
    class Polynomial {

    private:

        std::array<comp, N + 1> coeffs;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > N), int>::type = 0>
        Polynomial(std::array<comp, M> coeffs) noexcept;

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type = 0>
        Polynomial(Q... args) noexcept;

        auto value(comp x);

        auto solve();

        // TODO: Polynomial arithmetic
    };
}

#include <m/polynomial_impl.h>

#endif
