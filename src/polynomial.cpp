
#include <m/polynomial.h>

m::ComplexSolutions::ComplexSolutions() noexcept
    :inf(false) {}

m::ComplexSolutions::ComplexSolutions(std::unordered_set<std::complex<double>> finiteSet) noexcept
    :solutionSet(finiteSet), inf(false) {}

m::ComplexSolutions m::ComplexSolutions::empty() noexcept {

    return ComplexSolutions();
}

m::ComplexSolutions m::ComplexSolutions::finite(std::unordered_set<std::complex<double>> finiteSet) noexcept {

    return ComplexSolutions(finiteSet);
}

m::ComplexSolutions m::ComplexSolutions::infinite() noexcept {

    ComplexSolutions result;

    result.inf = true;

    return result;
}
