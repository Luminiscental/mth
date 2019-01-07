
#include <cmath>

#include <m/comp.h>

template <typename T>
bool m::util::checkZero(const T &x) {

    auto magnitude = std::abs(x);
    return magnitude <= m::EPSILON<decltype(magnitude)>;
}

template <typename T>
bool m::util::checkEqual(const T &a, const T &b) {

    return checkZero(a - b);
}
