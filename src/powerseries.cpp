
#include <mth/mth.h>

#include <mth/powerseries.h>

#include <mth/numeric.h>

mth::PowerSeries::PowerSeries() {

    generatingFunction = [] (size_t n) {

        return comp(0);
    };
}

mth::PowerSeries::PowerSeries(std::function<comp(size_t)> generatingFunction)
    :generatingFunction(generatingFunction) {}

mth::PowerSeries mth::PowerSeries::finite(const mth::Polynomial &equivalent) {

    const auto degree = equivalent.getDegree();

    if (degree.isInfinite()) return PowerSeries();

    const auto degreeValue = degree.getValue();

    auto gen = [&] (size_t n) {

        if (n > degreeValue) return comp(0);

        return equivalent.getCoeff(n);
    };

    return PowerSeries(gen);
}

mth::PowerSeries mth::PowerSeries::recursive(std::function<comp(comp)> recursion, const comp &constant) {

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

mth::comp mth::PowerSeries::getCoeff(size_t index) const {

    return generatingFunction(index);
}

mth::comp mth::PowerSeries::getPartial(const mth::comp &z, size_t index) const {

    if (util::isZero(z)) return getCoeff(0);

    mth::comp result;

    // Calculate from scratch
    if (lastPartialIndex == 0 || lastPartialIndex > index) {

        mth::comp term(1);

        for (size_t i = 0; i <= index; i++) {

            result += getCoeff(i) * term;
            term *= z;
        }

    // Add the extra terms
    } else if (lastPartialIndex < index) {

        using std::pow;

        result = lastPartial;
        mth::comp term = pow(z, lastPartialIndex + 1);

        for (size_t i = lastPartialIndex + 1; i <= index; i++) {

            result += getCoeff(i) * term;
            term *= z;
        }
    }

    lastPartialIndex = index;
    lastPartial = result;

    return result;
}

mth::comp mth::PowerSeries::getValue(const mth::comp &z) const {

    if (util::isZero(z)) return getCoeff(0);

    auto sequence = [&] (size_t n) {

        using std::pow;

        return getCoeff(n) * pow(z, n);
    };

    auto partialSequence = [&] (size_t n) {

        return getPartial(z, n);
    };

    return seriesLimit(partialSequence, sequence);
}

mth::PowerSeries mth::differentiate(const mth::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        auto c = series.getCoeff(index + 1);

        return c * comp(index + 1);
    };

    return PowerSeries(generatingFunction);
}

mth::PowerSeries mth::integrate(const mth::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        if (index == 0) return comp(0);

        auto c = series.getCoeff(index);

        return c / comp(index);
    };

    return PowerSeries(generatingFunction);
}

