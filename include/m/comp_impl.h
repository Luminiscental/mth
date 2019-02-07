
#include <m/m.h>
#include <m/vec.h>

template <>
bool m::util::isZero(const m::comp &x) {

    auto magnitude = std::abs(x);
    return magnitude <= m::EPSILON<decltype(magnitude)>;
}

#define BINDING(name, value)    template <typename T> const T &m::tcomp<T>:: name () const noexcept { return value ; } \
                                template <typename T>       T &m::tcomp<T>:: name ()       noexcept { return value ; }

BINDING(real, a)
BINDING(imag, b)

#undef BINDING

template <typename T>
m::tvec<T, 2> m::tcomp<T>::asCartesian() const noexcept {

    return tvec<T, 2>(a, b);
}

template <typename T>
m::tvec<T, 2> m::tcomp<T>::asPolar() const noexcept {

    return tvec<T, 2>(static_cast<T>(abs()), static_cast<T>(arg()));
}

template <typename T>
T m::tcomp<T>::absSqr() const noexcept {

    return a * a + b * b;
}

template <typename T>
double m::tcomp<T>::abs() const noexcept {

    auto ls = static_cast<double>(absSqr());
    
    if (util::isZero(ls)) return 0.0;

    return std::sqrt(ls);
}

template <typename T>
double m::tcomp<T>::arg() const noexcept {

    return std::atan2(static_cast<double>(b), static_cast<double>(a));
}

template <typename T>
m::tcomp<T> m::tcomp<T>::unit() const {

    auto l = abs();

    return *this / l;
}

template <typename T>
m::tcomp<T> m::tcomp<T>::conjugate() const noexcept {

    tcomp<T> result = a;

    result.b = -b;

    return result;
}

template <typename T>
m::tcomp<T> m::tcomp<T>::inverse() const {

    auto ls = absSqr();
    return conjugate() / ls;
}

template <typename T>
m::tcomp<T> m::tcomp<T>::rotation(const T &angle) {

    auto a = static_cast<double>(angle);

    auto c = static_cast<T>(std::cos(a));
    auto s = static_cast<T>(std::sin(a));

    return tcomp<T>(c, s);
}

template <typename T>
m::tcomp<T> m::tcomp<T>::fromCartesian(const tvec<T, 2> &vec) {

    return fromCartesian(vec.x(), vec.y());
}

template <typename T>
m::tcomp<T> m::tcomp<T>::fromPolar(const T &radius, const T &angle) {

    return radius * rotation(angle);
}

template <typename T>
m::tcomp<T> m::tcomp<T>::fromPolar(const m::tvec<T, 2> &polar) {

    return fromPolar(polar.x(), polar.y());
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator+=(const m::tcomp<T> &rhs) {

    a += rhs.a;
    b += rhs.b;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator+=(const T &rhs) {

    a += rhs;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator-=(const m::tcomp<T> &rhs) {

    a -= rhs.a;
    b -= rhs.b;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator-=(const T &rhs) {

    a -= rhs;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator*=(const m::tcomp <T> &rhs) {

    *this = *this * rhs;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator*=(const T &rhs) {

    a *= rhs;
    b *= rhs;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator/=(const m::tcomp<T> &rhs) {

    *this = *this / rhs;

    return *this;
}

template <typename T>
m::tcomp<T> &m::tcomp<T>::operator/=(const T &rhs) {

    a /= rhs;
    b /= rhs;

    return *this;
}

template <typename T>
m::tcomp<T> m::operator+(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result += rhs;
}

template <typename T>
m::tcomp<T> m::operator+(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result += rhs;
}

template <typename T>
m::tcomp<T> m::operator+(const T &lhs, const m::tcomp<T> &rhs) {

    return rhs + lhs;
}

template <typename T>
m::tcomp<T> m::operator-(const m::tcomp<T> &rhs) {

    tcomp<T> result = rhs;

    result.a = -result.a;
    result.b = -result.b;

    return result;
}

template <typename T>
m::tcomp<T> m::operator-(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tcomp<T> m::operator-(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tcomp<T> m::operator-(const T &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tcomp<T> m::operator*(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result;

    result.a = lhs.a * rhs.a - lhs.b * rhs.b;
    result.b = lhs.a * rhs.b + lhs.b * rhs.a;

    return result;
}

template <typename T>
m::tcomp<T> m::operator*(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result *= rhs;
}

template <typename T>
m::tcomp<T> m::operator*(const T &lhs, const m::tcomp<T> &rhs) {

    return rhs * lhs;
}

template <typename T>
m::tcomp<T> m::operator/(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
m::tcomp<T> m::operator/(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result /= rhs;
}

template <typename T>
m::tcomp<T> m::operator/(const T &lhs, const m::tcomp<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
bool m::operator==(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    return util::isEqual(lhs.real(), rhs.real()) && util::isEqual(lhs.imag(), rhs.imag());
}

template <typename T>
bool m::operator==(const m::tcomp<T> &lhs, const T &rhs) {

    return util::isZero(lhs.imag()) && util::isEqual(lhs.real(), rhs);
}

template <typename T>
bool m::operator==(const T &lhs, const m::tcomp<T> &rhs) {

    return rhs == lhs;
}

template <typename T>
bool m::operator!=(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    return !(lhs == rhs);
}

template <typename T>
bool m::operator!=(const m::tcomp<T> &lhs, const T &rhs) {

    return !(lhs == rhs);
}

template <typename T>
bool m::operator!=(const T &lhs, const m::tcomp<T> &rhs) {

    return !(lhs == rhs);
}

template <typename T>
std::ostream &m::operator<<(std::ostream &lhs, const m::tcomp<T> &rhs) {

    bool realZero = util::isZero(rhs.real());
    bool imagZero = util::isZero(rhs.imag());

    bool zero = realZero && imagZero;
    bool noneZero = !realZero && !imagZero;

    if (zero) return lhs << "0";

    // Bracket if there are multiple terms
    if (noneZero) lhs << "(";

    if (!realZero) lhs << rhs.real();

    if (!imagZero) {

        if (!realZero) {

            lhs << (rhs.imag() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.imag());

        } else {

            lhs << rhs.imag();
        }

        lhs << "i";
    }

    if (noneZero) lhs << ")";

    return lhs;
}

template <typename T>
double std::abs(const m::tcomp<T> &z) {

    return z.abs();
}

template <typename T>
m::tcomp<T> std::sqrt(const m::tcomp<T> &z) {

    auto p = z.asPolar();

    p.x() = sqrt(p.x());
    p.y() /= 2;

    return m::tcomp<T>::fromPolar(p);
}

template <typename T>
m::tcomp<T> std::exp(const m::tcomp<T> &z) {

    auto c = cos(z.imag());
    auto s = sin(z.imag());

    return exp(z.real()) * (c + m::i<T> * s);
}

template <typename T>
m::tcomp<T> std::log(const m::tcomp<T> &z) {

    auto p = z.asPolar();

    return log(p.x()) + m::tcomp<T>::fromCartesian(0, p.y());
}

template <typename T>
m::tcomp<T> std::cos(const m::tcomp<T> &z) {

    auto exponent = m::i<T> * z;

    return (exp(exponent) + exp(-exponent)) / static_cast<T>(2);
}

template <typename T>
m::tcomp<T> std::sin(const m::tcomp<T> &z) {

    auto exponent = m::i<T> * z;

    return (exp(exponent) - exp(-exponent)) / (2.0 * m::i<T>);
}

template <typename T>
m::tcomp<T> std::pow(const m::tcomp<T> &z, const m::tcomp<T> &exponent) {

    return exp(exponent * std::log(z));
}

template <typename T>
m::tcomp<T> std::pow(const m::tcomp<T> &z, const T &exponent) {

    return exp(exponent * std::log(exponent));
}

template <typename T>
m::tcomp<T> std::pow(const T &base, const m::tcomp<T> &z) {

    return exp(z * std::log(base));
}

template <typename T>
m::tcomp<T> std::pow(const m::tcomp<T> &z, size_t exponent) {

    auto result = z;

    for (size_t i = 1; i < exponent; i++) {

        result *= z;
    }

    return result;
}

template <typename T>
size_t std::hash<m::tcomp<T>>::operator()(const m::tcomp<T> &z) const {

    auto r = z.real();
    auto i = z.imag();

    return hash<decltype(r)>()(r) ^ hash<decltype(i)>()(i);
}
