
#include <mth/mth.h>

#include <mth/series.h>

#include <mth/numeric.h>

mth::Series::Series() {

    terms = [] (size_t index) {

        return (comp) 0;
    };
}

mth::Series::Series(std::function<comp(size_t)> terms)
    :terms(terms) {}

mth::comp mth::Series::getTerm(size_t index) const {

    return terms(index);
}

mth::comp mth::Series::getPartial(size_t index) const {

    mth::comp result;

    // Calculate from scratch
    if (lastPartialIndex == 0 || lastPartialIndex > index) {

        for (size_t i = 0; i <= index; i++) {

            result += getTerm(i);
        }
    } else if (lastPartialIndex < index) {

        result = lastPartial;
        
        for (size_t i = lastPartialIndex + 1; i <= index; i++) {

            result += getTerm(i);
        }
    }

    lastPartialIndex = index;
    lastPartial = result;

    return result;
}

mth::comp mth::Series::getLimit() const {

    auto partialSequence = [&] (size_t index) {

        return getPartial(index);
    };

    return seriesLimit(partialSequence, terms);
}

mth::Series mth::Series::recursive(std::function<mth::comp(mth::comp)> recursion, const mth::comp &init) {

    auto terms = [&] (size_t index) {

        comp result;

        for (int i = 0; i < index; i++) {

            result = recursion(result);
        }

        return result;
    };

    return Series(terms);
}

