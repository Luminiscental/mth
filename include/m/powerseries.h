#ifndef __m_series_h__
#define __m_series_h__

/* <m/series.h> - series header
 *      Defines the PowerSeries class, representing a complex power series.
 *      This has methods for evaluating at points, finding partial sums
 *      and full sums. It can be created either from a generating function,
 *      a finite polynomial, or a recursive coefficient relation.
 */

#include <functional>

#include <m/m.h>
#include <m/comp.h>
#include <m/polynomial.h>

namespace m {

    class PowerSeries {

    private:

        std::function<comp(size_t)> generatingFunction;

        // Lazy partial sum calculation
        mutable size_t lastPartialIndex;
        mutable comp lastPartial;

    public:

        // Default initializes to zero
        PowerSeries();
        PowerSeries(std::function<comp(size_t)> generatingFunction);

        comp getCoeff(size_t index) const;

        // TODO: Transforms (i.e. x -> (x - a))
        // TODO: Taylor series

        // Returns the partial sum up to index inclusive
        comp getPartial(const m::comp &z, size_t index) const;

        // Returns the numeric limit of the partial sums at z
        comp getValue(const m::comp &z) const;

        static PowerSeries finite(const Polynomial &equivalent);
        // TODO: Include n as a parameter
        static PowerSeries recursive(std::function<comp(comp)> recursion, const comp &constant);
    };

    PowerSeries differentiate(const PowerSeries &series);
    PowerSeries integrate(const PowerSeries &series);
}

#endif
