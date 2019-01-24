
#include <m/m.h>

#include <m/powerseries.h>

#include <m/numeric.h>

m::PowerSeries::PowerSeries() {

    generatingFunction = [] (size_t n) {

        return comp(0);
    };
}

m::PowerSeries::PowerSeries(std::function<comp(size_t)> generatingFunction)
    :generatingFunction(generatingFunction) {}

m::PowerSeries m::PowerSeries::finite(const m::Polynomial &equivalent) {

    const auto degree = equivalent.getDegree();

    if (degree.isInfinite()) return PowerSeries();

    const auto degreeValue = degree.getValue();

    auto gen = [&] (size_t n) {

        if (n > degreeValue) return comp(0);

        return equivalent.getCoeff(n);
    };

    return PowerSeries(gen);
}

m::PowerSeries m::PowerSeries::recursive(std::function<comp(comp)> recursion, const comp &constant) {

    auto generatingFunction = [&] (size_t index) {

        comp accumulate = constant;

        while (index > 0) {

            accumulate = recursion(accumulate);
            index--;
        }

        return accumulate;
    };

    return PowerSeries(generatingFunction);
}

m::comp m::PowerSeries::getCoeff(size_t index) const {

    return generatingFunction(index);
}

m::comp m::PowerSeries::getPartial(const m::comp &z, size_t index) const {

    if (util::isZero(z)) return getCoeff(0);

    m::comp result;

    // Calculate from scratch
    if (lastPartialIndex == 0 || lastPartialIndex > index) {

        m::comp term(1);

        for (size_t i = 0; i <= index; i++) {

            result += getCoeff(i) * term;
            term *= z;
        }

    // Add the extra terms
    } else if (lastPartialIndex < index) {

        result = lastPartial;
        m::comp term = std::pow(z, lastPartialIndex + 1);

        for (size_t i = lastPartialIndex + 1; i <= index; i++) {

            result += getCoeff(i) * term;
            term *= z;
        }
    }

    lastPartialIndex = index;
    lastPartial = result;

    return result;
}

m::comp m::PowerSeries::getValue(const m::comp &z) const {

    if (util::isZero(z)) return getCoeff(0);

    auto sequence = [&] (size_t n) {

        return getCoeff(n) * std::pow(z, n);
    };

    auto partialSequence = [&] (size_t n) {

        return getPartial(z, n);
    };

    return seriesLimit(partialSequence, sequence);
}

m::PowerSeries m::differentiate(const m::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        auto c = series.getCoeff(index + 1);

        return c * comp(index + 1);
    };

    return PowerSeries(generatingFunction);
}

m::PowerSeries m::integrate(const m::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        if (index == 0) return comp(0);

        auto c = series.getCoeff(index);

        return c / comp(index);
    };

    return PowerSeries(generatingFunction);
}

