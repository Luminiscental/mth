#ifndef __m_polynomial_h__
#define __m_polynomial_h__

#include <unordered_set>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#ifndef __m_constants_h__
#include <m/constants.h>
#endif

namespace std {

    template<>
    struct hash<std::complex<double>> {

        size_t operator()(const std::complex<double> &x) const {

            return std::hash<double>()(x.real()) ^ std::hash<double>()(x.imag());
        }
    };
}

namespace m {

    // NOTE: Nth degree polynomial
    
    class ComplexSolutions {

    private:

        std::unordered_set<std::complex<double>> solutionSet;
        bool inf;

        ComplexSolutions() noexcept;
        ComplexSolutions(std::unordered_set<std::complex<double>> finiteSet) noexcept;

    public:

        static ComplexSolutions empty() noexcept;

        static ComplexSolutions finite(std::unordered_set<std::complex<double>> finiteSet) noexcept;
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) noexcept {

            return ComplexSolutions(std::unordered_set<std::complex<double>>({args...}));
        }

        static ComplexSolutions infinite() noexcept;

        friend auto &operator<<(std::ostream &stream, const ComplexSolutions &solutions) {

            stream << std::fixed << std::setprecision(m_PRECISION);

            if (solutions.inf) {

                stream << "z in C";

            } else if (solutions.solutionSet.empty()) {

                stream << "no such z";

            } else {

                size_t n = solutions.solutionSet.size();
                size_t i = 0;

                for (auto root : solutions.solutionSet) {

                    if (i++ == n - 1) {

                        stream << "z = " << root;

                    } else {

                        stream << "z = " << root << ", or ";
                    }
                }
            }

            return stream;
        }
    };

    // TODO: Polynomial arithmetic

    template <size_t N>
    class Polynomial;

    template<>
    class Polynomial<0> {

    private:

        std::complex<double> coeff;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > 0), int>::type = 0>
        Polynomial(std::array<std::complex<double>, M> coeffs) noexcept
            :rootsValid(false), roots(ComplexSolutions::empty()) {
            
            coeff = coeffs[0];
        }

        Polynomial(std::complex<double> value) noexcept
            :coeff(value), rootsValid(false), roots(ComplexSolutions::empty()) {}

        auto value(std::complex<double> x) {

            return coeff;
        }

        auto solve() {

            if (rootsValid) return roots;

            rootsValid = true;

            if (util::checkZero(coeff)) return roots = ComplexSolutions::infinite();

            return roots = ComplexSolutions::empty();
        }

        friend auto &operator<<(std::ostream &stream, const Polynomial<0> &polynomial) {

            return stream << std::fixed << std::setprecision(m_PRECISION) << polynomial.coeff;
        }
    };

    template <size_t N>
    class Polynomial {

    private:

        std::array<std::complex<double>, N + 1> coeffs;
        ComplexSolutions roots;
        bool rootsValid;

    public:

        template <size_t M, typename std::enable_if<(M > N), int>::type = 0>
        Polynomial(std::array<std::complex<double>, M> coeffs)
            :rootsValid(false), roots(ComplexSolutions::empty()) {

            for (size_t i = 0; i < N + 1; i++) this->coeffs[i] = coeffs[i];
        }

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type = 0>
        Polynomial(Q... args)
            :Polynomial(std::array<std::complex<double>, N + 1> {args...}) {}

        auto value(std::complex<double> x) {

            std::complex<double> result = 0;
            std::complex<double> v = 1;

            for (size_t i = 0; i < N + 1; i++) {

                result += v * coeffs[i];
                v *= x;
            }

            return result;
        }

        auto solve() {

            if (rootsValid) return roots;

            rootsValid = true;

            if (util::checkZero(coeffs[N])) {

                return roots = Polynomial<N - 1>(coeffs).solve();
            }

            switch (N) {

                case 1: {

                    return roots = ComplexSolutions::finite(-coeffs[0] / coeffs[1]);
                }

                case 2: {

                        auto descriminant = coeffs[1] * coeffs[1] - 4.0 * coeffs[2] * coeffs[0];
                        auto offset = std::sqrt(descriminant);

                        auto lesser = -coeffs[1] - offset;
                        auto greater = -coeffs[1] + offset;

                        auto denom = 2.0 * coeffs[2];

                        if (util::checkEqual(lesser, greater)) return roots = ComplexSolutions::finite(lesser / denom);

                        return roots = ComplexSolutions::finite(lesser / denom, greater / denom);
                }

                default: {

                    // TODO: Numerical solution for higher degrees
                    // TODO: Cubics and quartics
                    
                    return roots = ComplexSolutions::empty();
                }
            }
        }

        friend auto &operator<<(std::ostream &lhs, const Polynomial<N> &rhs) {

            lhs << std::fixed << std::setprecision(m_PRECISION) << rhs.coeffs[0];

            if (N == 0) return lhs;
            
            lhs << " + " << rhs.coeffs[1] << "z";

            for (size_t i = 2; i < N + 1; i++) {

                lhs << " + " << rhs.coeffs[i];
                lhs << "z^" << i;
            }

            return lhs;
        }

    };
}

#endif
