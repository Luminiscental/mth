#ifndef mth_powerseries_h__
#define mth_powerseries_h__

/* <mth/powerseries.h> - power series header
 *      Defines the PowerSeries class, representing a complex power series.
 *      This has methods for evaluating at points, finding partial sums
 *      and full sums. It can be created either from a generating function,
 *      a finite polynomial, or a recursive coefficient relation.
 */

#include <functional>

#include <mth/mth.h>
#include <mth/comp.h>
#include <mth/polynomial.h>
#include <mth/series.h>

namespace mth {

    class PowerSeries {

    private:

            std::function<comp(size_t)> generatingFunction;

    public:

        // Default initializes to zero
        PowerSeries();
        PowerSeries(std::function<comp(size_t)> generatingFunction);

        comp getCoeff(size_t index) const;

        // TODO: Transforms (i.e. x -> (x - a))
        // TODO: Taylor series
        
        // Returns the series evaluated at z
        Series series(const comp &z) const;

        static PowerSeries finite(const Polynomial &equivalent);

        // TODO: Include n as a parameter
        static PowerSeries recursive(std::function<comp(comp)> recursion, const comp &constant);
    };

    PowerSeries differentiate(const PowerSeries &series);
    PowerSeries integrate(const PowerSeries &series);
}

#endif
