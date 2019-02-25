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

            std::function<comp(size_t)> generatingFunction = [] (size_t n) {

                return comp(0);
            };

            bool isTrivial = false;
            Polynomial trivialSeries;

    public:

        // Default initializes to zero
        PowerSeries() = default;

        // Initialize from a generating function
        PowerSeries(std::function<comp(size_t)> generatingFunction);

        // Get the coefficient at an index
        comp getCoeff(size_t index) const;

        // TODO: Transforms (i.e. x -> (x - a))
        // TODO: Taylor series
        
        // Returns the series evaluated at z
        Series series(const comp &z) const;

        // Create a finite power series from a polynomial
        // TODO: Have this as a cast constructor too maybe
        static PowerSeries finite(const Polynomial &equivalent);

        // Create a power series from a recursive relation of the coefficients
        // TODO: Include n as a parameter
        static PowerSeries recursive(std::function<comp(comp)> recursion, const comp &constant);
    };

    // Differentiation and integration of power series
    
    PowerSeries differentiate(const PowerSeries &series);
    PowerSeries integrate(const PowerSeries &series);
}

#endif
