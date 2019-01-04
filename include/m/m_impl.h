#ifdef __m_impl__
#undef __m_impl__

#include <cmath>
#include <limits>

template <typename T>
constexpr T m::EPSILON = std::numeric_limits<T>::epsilon();

template <typename T>
constexpr T m::PI = static_cast<T>(3.14159265358979323846264338327950288419716939937510582097494459230781640629);
 
template <typename T>
constexpr T m::TAU = static_cast<T>(6.28318530717958647692528676655900576839433879875021164194988918461563281257);

// NOTE: Forward def for use of std::abs

namespace m {

    template <typename T>
    class tcomp;
}

namespace std {

    template <typename T>
    double abs(const m::tcomp<T> &z);
}

template <typename T>
auto m::util::checkZero(const T &x) {

    using std::abs;

    auto magnitude = abs(x);
    return magnitude <= m::EPSILON<decltype(magnitude)>;
}

template <typename T>
auto m::util::checkEqual(const T &a, const T &b) {

    return checkZero(a - b);
}

#define __m_impl__
#endif
