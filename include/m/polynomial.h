#ifndef __m_polynomial_h__
#define __m_polynomial_h__

#include <unordered_set>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>

#ifndef __m_complex_h__
#include <m/complex.h>
#endif

// TODO: Move to complex and add for other classes
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
    auto &operator<<(std::ostream &lhs, const Polynomial<N> &rhs);

    // NOTE: Nth degree polynomial
    
    class ComplexSolutions {

    private:

        std::unordered_set<m::comp> solutionSet;
        bool inf;

        ComplexSolutions() noexcept;
        ComplexSolutions(std::unordered_set<m::comp> finiteSet) noexcept;

    public:

        static ComplexSolutions empty() noexcept;

        static ComplexSolutions finite(std::unordered_set<m::comp> finiteSet) noexcept;
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept;

        static ComplexSolutions infinite() noexcept;
    };

    // TODO: Polynomial arithmetic

    template <size_t N>
    class Polynomial;

    template<>
    class Polynomial<0> {

    private:

        m::comp coeff;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > 0), int>::type = 0>
        Polynomial(std::array<m::comp, M> coeffs) noexcept;

        Polynomial(m::comp value) noexcept;

        auto value(m::comp x);

        auto solve();
    };

    template <size_t N>
    class Polynomial {

    private:

        std::array<m::comp, N + 1> coeffs;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > N), int>::type = 0>
        Polynomial(std::array<m::comp, M> coeffs);

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type = 0>
        Polynomial(Q... args);

        auto value(m::comp x);

        auto solve();
    };
}

#include <m/polynomial_impl.h>

#endif
