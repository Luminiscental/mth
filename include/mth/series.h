#ifndef mth_series_h__
#define mth_series_h__

#include <functional>

#include <mth/comp.h>

namespace mth {

    class Series {

    private:

        std::function<comp(size_t)> terms;

        mutable size_t lastPartialIndex;
        mutable comp lastPartial;

    public:

        // Default initializes to zero
        Series();
        Series(std::function<comp(size_t)> terms);

        comp getTerm(size_t index) const;

        // Returns the partial sum up to index inclusive
        comp getPartial(size_t index) const;

        // Returns the numeric limit of the partial sums at z
        comp getLimit() const;

        template <typename ...Q>
        static Series finite(Q... terms) {

            std::vector<comp> nonTrivialTerms {static_cast<comp>(terms) ...};

            return Series([nonTrivialTerms] (size_t index) {

                if (index < nonTrivialTerms.size()) {

                    return nonTrivialTerms[index];

                } else {

                    return (comp) 0;
                }
            });
        }

        static Series recursive(std::function<comp(comp)> recursion, const comp &init);
    };
}

#endif
