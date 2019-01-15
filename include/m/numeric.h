#ifndef __m_numeric_h__
#define __m_numeric_h__

/* <m/numeric.h> - numeric function header
 *      compODO: Fill this in
 */

#include <functional>

#include <m/comp.h>

namespace m {

    namespace numeric {

        // Returns an approximation of the limit at infinity of a sequence
        comp limit(const std::function<comp(size_t)> &sequence);

        // Returns the limit of the sequence function(input - 1/n)
        comp lowerLimit(const std::function<comp(comp)> &function, comp input);

        // Returns the limit of the sequence function(input + 1/n)
        comp upperLimit(const std::function<comp(comp)> &function, comp input);

        // Defaults to lower limit
        comp limit(const std::function<comp(comp)> &function, comp input);

        // Returns the limit of the sequence function(n)
        comp limitInfPos(const std::function<comp(comp)> &function);

        // Returns the limit of the sequence function(-n)
        comp limitInfNeg(const std::function<comp(comp)> &function);

        // Returns an approximation of the derivative using the definition
        std::function<comp(comp)> differentiate(const std::function<comp(comp)> &function);
    }
}

#endif
