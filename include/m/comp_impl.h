#ifdef __m_impl__
#undef __m_impl__

#define BINDING(name, value)   template <typename T> const auto &m::tcomp<T>:: name () const noexcept { return value ; } \
                        template <typename T> auto &m::tcomp<T>:: name () noexcept { return value ; }

BINDING(real, a)
BINDING(imag, b)

#undef BINDING

template <typename T>
m::tvec<T, 2> m::tcomp<T>::asPolar() const noexcept {

    return tvec<T, 2>(static_cast<T>(abs()), static_cast<T>(arg()));
}

template <typename T>
m::tcomp<T>::operator m::tvec<T, 2>() noexcept {

    return tvec<T, 2>(a, b);
}

template <typename T>
auto m::tcomp<T>::absSqr() const noexcept {

    return a * a + b * b;
}

template <typename T>
auto m::tcomp<T>::abs() const noexcept {

    auto ls = static_cast<double>(absSqr());
    
    if (util::checkZero(ls)) return 0.0;

    return std::sqrt(ls);
}

template <typename T>
auto m::tcomp<T>::arg() const noexcept {

    return std::atan2(static_cast<double>(a), static_cast<double>(b));
}

template <typename T>
auto m::tcomp<T>::unit() const {

    auto l = abs();

    if (util::checkZero(l)) throw std::invalid_argument("m::exception: tcomp zero has no unit equivalent");

    return *this / l;
}

template <typename T>
auto m::tcomp<T>::conjugate() const noexcept {

    tcomp<T> result = a;

    result.b = -b;

    return result;
}

template <typename T>
auto m::tcomp<T>::inverse() const {

    auto ls = absSqr();

    if (util::checkZero(ls)) throw std::invalid_argument("m::exception: tcomp zero has no inverse");

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
m::tcomp<T> m::tcomp<T>::fromPolar(const T &radius, const T &angle) {

    return radius * rotation(angle);
}

template <typename T>
m::tcomp<T> m::tcomp<T>::fromPolar(const m::tvec<T, 2> &polar) {

    return fromPolar(polar.x(), polar.y());
}

template <typename T>
auto &m::tcomp<T>::operator+=(const m::tcomp<T> &rhs) {

    a += rhs.a;
    b += rhs.b;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator+=(const T &rhs) {

    a += rhs;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator-=(const m::tcomp<T> &rhs) {

    a -= rhs.a;
    b -= rhs.b;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator-=(const T &rhs) {

    a -= rhs;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator*=(const m::tcomp <T> &rhs) {

    *this = *this * rhs;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator*=(const T &rhs) {

    a *= rhs;
    b *= rhs;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator/=(const m::tcomp<T> &rhs) {

    *this = *this / rhs;

    return *this;
}

template <typename T>
auto &m::tcomp<T>::operator/=(const T &rhs) {

    a /= rhs;
    b /= rhs;

    return *this;
}

template <typename T>
auto m::operator+(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result += rhs;
}

template <typename T>
auto m::operator+(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result += rhs;
}

template <typename T>
auto m::operator+(const T &lhs, const m::tcomp<T> &rhs) {

    return rhs + lhs;
}

template <typename T>
auto m::operator-(const m::tcomp<T> &rhs) {

    tcomp<T> result = rhs;

    result.a = -result.a;
    result.b = -result.b;

    return result;
}

template <typename T>
auto m::operator-(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator-(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator-(const T &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator*(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    tcomp<T> result;

    result.a = lhs.a * rhs.a - lhs.b * rhs.b;
    result.b = lhs.a * rhs.b + lhs.b * rhs.a;

    return result;
}

template <typename T>
auto m::operator*(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result *= rhs;
}

template <typename T>
auto m::operator*(const T &lhs, const m::tcomp<T> &rhs) {

    return rhs * lhs;
}

template <typename T>
auto m::operator/(const m::tcomp<T> &lhs, const m::tcomp<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
auto m::operator/(const m::tcomp<T> &lhs, const T &rhs) {

    tcomp<T> result = lhs;

    return result /= rhs;
}

template <typename T>
auto m::operator/(const T &lhs, const m::tcomp<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
auto &m::operator<<(std::ostream &lhs, const m::tcomp<T> &rhs) {

    bool nonZero = false;

    lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

    if (!util::checkZero(rhs.real())) {

        lhs << rhs.real();
        nonZero = true;
    }

    if (!util::checkZero(rhs.imag())) {

        if (nonZero) {

            lhs << (rhs.imag() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.imag());

        } else {

            lhs << rhs.imag();
            nonZero = true;
        }

        lhs << "i";
    }

    if (!nonZero) lhs << "0";

    return lhs << ")";
}

auto m::operator "" _i(long double im) {

    return m::tcomp<float>(static_cast<float>(im));
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

    auto c = static_cast<T>(std::cos(z.imag()));
    auto s = static_cast<T>(std::sin(z.imag()));

    return exp(z.real()) * (c + m::i * s);
}

template <typename T>
m::tcomp<T> std::cos(const m::tcomp<T> &z) {

    return exp(m::i * z).real();
}

template <typename T>
m::tcomp<T> std::sin(const m::tcomp<T> &z) {

    return exp(m::i * z).imag();
}

#define __m_impl__
#endif
