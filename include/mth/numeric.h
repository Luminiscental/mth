#ifndef mth_numeric_h__
#define mth_numeric_h__

/* <mth/numeric.h> - numeric function header
 *      Defines functions to numerically estimate the value of limits
 *      and derivatives. Input functions can be any std::function on 
 *      mth::comp but non-well-behaved functions may cause issues.
 */

#include <functional>

#include <mth/mth.h>
#include <mth/comp.h>

// TODO: Limit of recursive sequence
// TODO: Less naive methods

namespace mth {

    // Returns an approximation of the limit at infinity of a sequence
    comp limit(const std::function<comp(size_t)> &sequence);

    // Returns an approximation of the limit at infinity of a series
    comp seriesLimit(const std::function<comp(size_t)> &partialSum, const std::function<comp(size_t)> &sequence);

    // Returns the limit of the sequence approaching from below (parallel with real axis)
    comp lowerLimit(const std::function<comp(comp)> &function, const comp &input);

    // Returns the limit of the sequence approaching from above (parallel with real axis)
    comp upperLimit(const std::function<comp(comp)> &function, const comp &input);

    // Defaults to lower limit
    comp limit(const std::function<comp(comp)> &function, const comp &input);

    // Returns the limit of the sequence function(n)
    comp limitInfPos(const std::function<comp(comp)> &function);

    // Returns the limit of the sequence function(-n)
    comp limitInfNeg(const std::function<comp(comp)> &function);

    // Returns an approximation of the derivative using numeric::limit
    std::function<comp(comp)> differentiate(const std::function<comp(comp)> &function);
}

#endif
