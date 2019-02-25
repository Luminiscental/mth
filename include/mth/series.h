#ifndef mth_series_h__
#define mth_series_h__

// TODO: Comment

#include <functional>

#include <mth/comp.h>

namespace mth {

    class Series {

    private:

        std::function<comp(size_t)> terms;

        mutable size_t lastPartialIndex;
        mutable comp lastPartial;

        bool isTrivial = false;
        comp trivialSum;

    public:

        // Default initializes to zero
        Series();
        Series(std::function<comp(size_t)> terms);

        comp getTerm(size_t index) const;

        // Returns the partial sum up to index inclusive
        comp getPartial(size_t index) const;

        // Returns the numeric limit of the partial sums at z
        comp getLimit() const;

        static Series finite(std::vector<comp> terms);

        template <typename ...Q>
        static Series finite(Q... terms) {

            std::vector<comp> nonTrivialTerms {static_cast<comp>(terms) ...};

            Series result([nonTrivialTerms] (size_t index) {

                if (index < nonTrivialTerms.size()) {

                    return nonTrivialTerms[index];

                } else {

                    return (comp) 0;
                }
            });

            result.isTrivial = true;
            result.trivialSum = result.getPartial(nonTrivialTerms.size());
            
            return result;
        }

        static Series recursive(std::function<comp(comp)> recursion, const comp &init);
    };
}

#endif
