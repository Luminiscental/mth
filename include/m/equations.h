#ifndef __m_equations_h__
#define __m_equations_h__

#include <set>
#include <cmath>
#include <complex>

#ifndef __m_vec_h__
#include <m/vec.h>
#endif

#ifndef __m_mat_h__
#include <m/mat.h>
#endif 

#ifndef __m_constants_h__
#include <m/constants.h>
#endif

namespace m {

    namespace eq {

        // TODO: Degenerate cases

        template <typename Solution>
        class System {

        public:

            virtual Solution solve() = 0;
        };

        template <typename T, size_t N>
        struct VectorComparator {

            bool operator() (const tvec<T, N> &a, const tvec<T, N> &b) {

                return (a.x() < b.x()) || (a.x() == b.x() && a.y() < b.y());
            }
        };
        
        // Linear system of equations in N variables (finding solution)
        
        template <typename T, size_t N>
        class LinearSystem : public System<std::set<tvec<T, N>, VectorComparator<T, N>>> {

        private:

            const tvec<T, N> aux;
            const tmat<T, N> coeffs;
            const bool singular;

        public:

            LinearSystem(tmat<T, N> coeffs, tvec<T, N> aux) : coeffs(coeffs), aux(aux), singular(util::checkZero(coeffs.det())) {}

            std::set<tvec<T, N>, VectorComparator<T, N>> solve() override {

                std::set<tvec<T, N>, VectorComparator<T, N>> result;

                if (!singular) {

                    result.insert(coeffs.inverse() * aux);

                } else {


                }

                return result;
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

                // TODO: Numerical solution for higher degrees
                // TODO: Cubics and quartics

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
                if (!util::checkEqual(lesser, greater)) result.insert(greater / denom);

                return result;
            }
        };
    }
}

#endif
