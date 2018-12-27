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

        ComplexSolutions() : inf(false) {}
        ComplexSolutions(std::unordered_set<std::complex<double>> finiteSet) : solutionSet(finiteSet), inf(false) {}

    public:

        static ComplexSolutions empty();

        static ComplexSolutions finite(std::unordered_set<std::complex<double>> finiteSet);
        
        template <typename ...Q>
        static ComplexSolutions finite(Q... args) {

            return ComplexSolutions(std::unordered_set<std::complex<double>>({args...}));
        }

        static ComplexSolutions infinite();

        friend std::ostream &operator<<(std::ostream &stream, const ComplexSolutions &solutions) {

            stream << std::fixed << std::setprecision(

#ifdef m_PRECISION

            m_PRECISION

#else

            2

#endif

            );

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

    template <size_t N>
    class Polynomial;

    template<>
    class Polynomial<0> {

    private:

        std::complex<double> coeff;

    public:

        template <size_t M, typename std::enable_if<(M > 0), int>::type = 0>
        Polynomial(std::array<std::complex<double>, M> coeffs) {
            
            coeff = coeffs[0];
        }

        Polynomial(std::complex<double> value) : coeff(value) {}

        std::complex<double> value(std::complex<double> x) {

            return coeff;
        }

        ComplexSolutions solve() {

            if (util::checkZero(coeff)) return ComplexSolutions::infinite();

            return ComplexSolutions::empty();
        }

        friend std::ostream &operator<<(std::ostream &stream, const Polynomial<0> &polynomial) {

            stream << std::fixed << std::setprecision(

#ifdef m_PRECISION

            m_PRECISION

#else

            2

#endif

            );

            return stream << polynomial.coeff;
        }
    };

    template <size_t N>
    class Polynomial {

    private:

        std::array<std::complex<double>, N + 1> coeffs;

    public:

        template <size_t M, typename std::enable_if<(M > N), int>::type = 0>
        Polynomial(std::array<std::complex<double>, M> coeffs) {

            for (size_t i = 0; i < N + 1; i++) this->coeffs[i] = coeffs[i];
        }

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type = 0>
        Polynomial(Q... args) : Polynomial(std::array<std::complex<double>, N + 1> {args...}) {}

        std::complex<double> value(std::complex<double> x) {

            std::complex<double> result = 0;
            std::complex<double> v = 1;

            for (size_t i = 0; i < N + 1; i++) {

                result += v * coeffs[i];
                v *= x;
            }

            return result;
        }

        ComplexSolutions solve() {

            if (util::checkZero(coeffs[N])) {

                return Polynomial<N - 1>(coeffs).solve();
            }

            switch (N) {

                case 1: {

                    return ComplexSolutions::finite(-coeffs[0] / coeffs[1]);
                }

                case 2: {

                        std::complex<double> descriminant = coeffs[1] * coeffs[1] - 4.0 * coeffs[2] * coeffs[0];
                        std::complex<double> offset = std::sqrt(descriminant);

                        std::complex<double> lesser = -coeffs[1] - offset;
                        std::complex<double> greater = -coeffs[1] + offset;

                        std::complex<double> denom = 2.0 * coeffs[2];

                        if (util::checkEqual(lesser, greater)) return ComplexSolutions::finite(lesser / denom);

                        return ComplexSolutions::finite(lesser / denom, greater / denom);
                }

                default: {

                    // TODO: Numerical solution for higher degrees
                    // TODO: Cubics and quartics
                    
                    return ComplexSolutions::empty();
                }
            }
        }

        friend std::ostream &operator<<(std::ostream &stream, const Polynomial<N> &polynomial) {

            stream << std::fixed << std::setprecision(

#ifdef m_PRECISION

            m_PRECISION

#else

            2

#endif

            );

            stream << polynomial.coeffs[0];
            if (N == 0) return stream;
            
            stream << " + " << polynomial.coeffs[1] << "z";

            for (size_t i = 2; i < N + 1; i++) {

                stream << " + " << polynomial.coeffs[i];
                stream << "z^" << i;
            }

            return stream;
        }

    };
}

#endif
