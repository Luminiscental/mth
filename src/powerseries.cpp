
#include <mth/mth.h>

#include <mth/powerseries.h>

#include <mth/numeric.h>

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

    PowerSeries result{gen};

    result.isTrivial = true;
    result.trivialSeries = equivalent;

    return result;
}

mth::PowerSeries mth::PowerSeries::recursive(std::function<comp(comp)> recursion, const comp &constant) {

    auto generatingFunction = [&] (size_t index) {

        auto accumulate = constant;

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

mth::Series mth::PowerSeries::series(const mth::comp &z) const {

    if (isTrivial) {

        std::vector<comp> result;

        auto deg = trivialSeries.getDegree();
        if (deg.isInfinite()) return mth::Series::finite(0);

        auto acc = comp{1};

        for (size_t i = 0; i <= deg.getValue(); i++) {

            auto coeff = trivialSeries.getCoeff(i);
            result.push_back(coeff * acc);

            acc *= z;
        }

        return mth::Series::finite(result);

    } else {

        return Series([&] (size_t index) {

            using std::pow;
            return generatingFunction(index) * pow(z, index);

        });
    }
}

mth::PowerSeries mth::differentiate(const mth::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        auto c = series.getCoeff(index + 1);

        return c * comp{static_cast<double>(index + 1)};
    };

    return PowerSeries(generatingFunction);
}

mth::PowerSeries mth::integrate(const mth::PowerSeries &series) {

    auto generatingFunction = [&] (size_t index) {

        if (index == 0) return comp{0};

        auto c = series.getCoeff(index);

        return c / comp{static_cast<double>(index)};
    };

    return PowerSeries(generatingFunction);
}

