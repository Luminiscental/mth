#ifndef __Maths_equations_h__
#define __Maths_equations_h__

#ifndef __Maths_vec_h__
#include <Maths/vec.h>
#endif

#ifndef __Maths_mat_h__
#include <Maths/mat.h>
#endif 

#ifndef __Maths_constants_h__
#include <Maths/constants.h>
#endif

#include <set>
#include <cmath>
#include <complex>

namespace m {

    namespace eq {

        template <typename Solution>
        class System {

        public:
            virtual Solution solve() = 0;
        };
        
        // Linear system of equations in N variables (finding solution)
        
        template <typename T, size_t N>
        class LinearSystem : public System<tvec<T, N>> {

        private:
            tvec<T, N> aux;
            tmat<T, N> coeffs;

        public:
            LinearSystem(tmat<T, N> coeffs, tvec<T, N> aux) : coeffs(coeffs), aux(aux) {}

            tvec<T, N> solve() override {

                return coeffs.inverse() * aux;
            }
        };

        struct ComplexComparator {

            bool operator() (const std::complex<double> &a, const std::complex<double> &b) {

                return (a.real() < b.real()) || (a.real() == b.real() && a.imag() < b.imag());
            }
        };

        // Nth degree polynomial over C (finding roots)

        template <size_t N>
        class Polynomial : public System<std::set<std::complex<double>, ComplexComparator>> {

        private:
            tvec<std::complex<double>, N + 1> coeffs;

        public:
            Polynomial(tvec<std::complex<double>, N + 1> coeffs) : coeffs(coeffs) {}

            std::complex<double> value(std::complex<double> x) {

                std::complex<double> result = 0;
                std::complex<double> v = 1;

                for (size_t i = 0; i < N + 1; i++) {

                    result += v * coeffs.get(i);
                    v *= x;
                }

                return result;
            }

            std::set<std::complex<double>, ComplexComparator> solve() override {

                std::set<std::complex<double>, ComplexComparator> result;

                return result;
            }
        };

        template<>
        class Polynomial<1> : public System<std::set<std::complex<double>, ComplexComparator>> {

        private:
            tvec<std::complex<double>, 2> coeffs;
            std::complex<double> root;
        
        public:
            Polynomial(tvec<std::complex<double>, 2> coeffs) : coeffs(coeffs) {

                root = -coeffs.get(0) / coeffs.get(1);
            }

            std::complex<double> value(std::complex<double> x) {

                return coeffs.get(0) + x * coeffs.get(1);
            }

            std::set<std::complex<double>, ComplexComparator> solve() override {

                std::set<std::complex<double>, ComplexComparator> result;

                result.insert(root);

                return result;
            }
        };

        template<>
        class Polynomial<2> : public System<std::set<std::complex<double>, ComplexComparator>> {

        private:
            tvec<std::complex<double>, 3> coeffs;
            std::complex<double> descriminant;

        public:
            Polynomial(tvec<std::complex<double>, 3> coeffs) : coeffs(coeffs) {

                descriminant = coeffs.get(1) * coeffs.get(1) - 4.0 * coeffs.get(2) * coeffs.get(0);
            }

            std::complex<double> value(std::complex<double> x) {

                return coeffs.get(0) + x * (coeffs.get(1) + x * coeffs.get(2));
            }

            std::set<std::complex<double>, ComplexComparator> solve() override {

                std::complex<double> offset = std::sqrt(descriminant);

                std::complex<double> lesser = -coeffs.get(1) - offset;
                std::complex<double> greater = -coeffs.get(1) + offset;

                std::set<std::complex<double>, ComplexComparator> result;

                std::complex<double> denom = 2.0 * coeffs.get(2);

                result.insert(lesser / denom);
                result.insert(greater / denom);

                return result;
            }
        };
    }
}

#endif
