
#include <m/polynomial.h>

m::ComplexSolutions m::ComplexSolutions::empty() {

    return ComplexSolutions();
}

m::ComplexSolutions m::ComplexSolutions::finite(std::unordered_set<std::complex<double>> finiteSet) {

    return ComplexSolutions(finiteSet);
}

m::ComplexSolutions m::ComplexSolutions::infinite() {

    ComplexSolutions result;

    result.inf = true;

    return result;
}
