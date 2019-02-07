
#include <cmath>

template <typename T>
bool m::util::isZero(const T &x) {

    auto magnitude = std::abs(x);
    return magnitude <= m::EPSILON<decltype(magnitude)>;
}

template <typename T>
bool m::util::isEqual(const T &a, const T &b) {

    return isZero(a - b);
}
