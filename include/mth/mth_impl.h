
#include <cmath>

template <typename T>
bool mth::util::isZero(const T &x) {

    using std::abs;
    auto magnitude = abs(x);
    return magnitude <= mth::EPSILON<decltype(magnitude)>;
}

template <typename T>
bool mth::util::isEqual(const T &a, const T &b) {

    return isZero(a - b);
}
