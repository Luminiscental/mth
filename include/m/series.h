#ifndef __m_series_h__
#define __m_series_h__

/* <m/series.h> - series header
 *      TODO: Fill out
 */

#include <functional>

#include <m/comp.h>
#include <m/polynomial.h>

namespace m {

    class Series {

    private:

        std::function<comp(size_t)> generatingFunction;

        // Lazy partial sum calculation
        mutable size_t lastPartialIndex;
        mutable comp lastPartial;

    public:

        // Default initializes to zero
        Series();
        Series(std::function<comp(size_t)> generatingFunction);

        comp getCoeff(size_t index) const;

        // TODO: Transforms (i.e. x -> (x - a))
        // TODO: Taylor series

        // Returns the partial sum up to index inclusive
        comp getPartial(const m::comp &z, size_t index) const;

        // Returns the numeric limit of the partial sums at z
        comp getValue(const m::comp &z) const;

        static Series finite(const Polynomial &equivalent);
        // TODO: Include n as a parameter
        static Series recursive(std::function<comp(comp)> recursion, const comp &constant);
    };

    Series differentiate(const Series &series);
    Series integrate(const Series &series);
}

#endif
