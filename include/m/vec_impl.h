#ifdef __m_impl__
#undef __m_impl__

#include <cmath>

template <typename T, size_t N>
m::tvec<T, N>::tvec(const std::array<T, N> &values) noexcept
    :values(values) {}

template <typename T, size_t N>
auto m::tvec<T, N>::begin() {

    return values.begin();
}

template <typename T, size_t N>
auto m::tvec<T, N>::begin() const {

    return cbegin();
}

template <typename T, size_t N>
auto m::tvec<T, N>::end() {

    return values.end();
}

template <typename T, size_t N>
auto m::tvec<T, N>::end() const {

    return cend();
}

template <typename T, size_t N>
auto m::tvec<T, N>::cbegin() const {

    return values.cbegin();
}

template <typename T, size_t N>
auto m::tvec<T, N>::cend() const {

    return values.cend();
}

template <typename T, size_t N>
auto m::tvec<T, N>::rbegin() {

    return values.rbegin();
}

template <typename T, size_t N>
auto m::tvec<T, N>::rend() {

    return values.rend();
}

template <typename T, size_t N>
auto m::tvec<T, N>::crbegin() const {

    return values.crbegin();
}

template <typename T, size_t N>
auto m::tvec<T, N>::crend() const {

    return values.crend();
}

template <typename T, size_t N>
T &m::tvec<T, N>::get(size_t index) {

    if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

    return values[index];
}

template <typename T, size_t N>
const T &m::tvec<T, N>::get(size_t index) const {

    if (index > N - 1) throw std::out_of_range("m::exception: vector accessed out of range");

    return values[index];
}

#define BINDING(name, value) template <typename T, size_t N> const T &m::tvec<T, N>:: name () const { return value ; } \
                                   template <typename T, size_t N> T &m::tvec<T, N>:: name ()       { return value; }

    BINDING(x, get(0))
    BINDING(y, get(1))
    BINDING(z, get(2))
    BINDING(w, get(3))

#undef BINDING

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 1), int>::type>
m::tvec<T, 2> m::tvec<T, N>::xy() const {

    return tvec<T, 2>(x(), y());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 1), int>::type>
m::tvec<T, 2> m::tvec<T, N>::yz() const {

    return tvec<T, 2>(y(), z());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 1), int>::type>
m::tvec<T, 2>  m::tvec<T, N>::zw() const {

    return tvec<T, 2>(z(), w());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 2), int>::type>
m::tvec<T, 3>  m::tvec<T, N>::xyz() const {

    return tvec<T, 3>(x(), y(), z());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 2), int>::type>
m::tvec<T, 3>  m::tvec<T, N>::yzw() const {

    return tvec<T, 3>(y(), z(), w());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 3), int>::type>
m::tvec<T, 4> m::tvec<T, N>::xyzw() const {

    return tvec<T, 4>(x(), y(), z(), w());
}

template <typename T, size_t N>
auto m::tvec<T, N>::magnSqr() const noexcept {

    auto result = static_cast<T>(0);

    for (const auto &value : *this) {

        result += value * value;
    }

    return result;
}

template <typename T, size_t N>
auto m::tvec<T, N>::magn() const noexcept {

    auto ls = static_cast<double>(magnSqr());

    if (util::checkZero(ls)) return 0.0;

    return std::sqrt(ls);
}

template <typename T, size_t N>
auto m::tvec<T, N>::dot(const m::tvec<T, N> &rhs) const {

    T result = 0;

    for (size_t i = 0; i < N; i++) {

        result += this->get(i) * rhs.get(i);
    }

    return result;
}

template <typename T, size_t N>
auto m::tvec<T, N>::unit() const {

    auto l = static_cast<T>(magn());

    if (util::checkZero(l)) throw std::invalid_argument("m::exception: unit() called on zero vector");

    return *this / l;
}

template <typename T, size_t N>
auto m::tvec<T, N>::dot(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) {

    return lhs.dot(rhs);
}

template <typename T, size_t N>
auto &m::tvec<T, N>::operator+=(const m::tvec<T, N> &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) += rhs.get(i);
    }

    return *this;
}

template <typename T, size_t N>
auto &m::tvec<T, N>::operator-=(const m::tvec<T, N> &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) -= rhs.get(i);
    }

    return *this;
}

template <typename T, size_t N>
auto &m::tvec<T, N>::operator*=(const T &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) *= rhs;
    }

    return *this;
}

template <typename T, size_t N>
auto &m::tvec<T, N>::operator/=(const T &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) /= rhs;
    }

    return *this;
}

template <typename T>
auto m::vec::cross(const m::tvec<T, 3> &lhs, const m::tvec<T, 3> &rhs) {

    return tvec<T, 3>(lhs.get(2) * rhs.get(3) - lhs.get(3) * rhs.get(2),
                      lhs.get(3) * rhs.get(1) - lhs.get(1) * rhs.get(3),
                      lhs.get(1) * rhs.get(2) - lhs.get(2) * rhs.get(1));
}

template <typename T, size_t N>
auto m::operator+(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) {

    auto result = lhs;

    return result += rhs;
}

template <typename T, size_t N>
auto m::operator-(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) {

    auto result = lhs;

    return result -= rhs;
}

template <typename T, size_t N>
auto m::operator-(const m::tvec<T, N> &rhs) {

    auto result = rhs;

    for (size_t i = 0; i < N; i++) {

        result.get(i) = -rhs.get(i);
    }

    return result;
}

template <typename T, size_t N>
auto m::operator*(const T &lhs, const m::tvec<T, N> &rhs) {

    auto result = rhs;

    return result *= lhs;
}

template <typename T, size_t N>
auto m::operator*(const m::tvec<T, N> &lhs, const T &rhs) {

    return rhs * lhs;
}

template <typename T, size_t N>
auto m::operator/(const m::tvec<T, N> &lhs, const T &rhs) {

    tvec<T, N> result = lhs;

    return result /= rhs;
}

template <typename T, size_t N>
auto &m::operator<<(std::ostream &lhs, const m::tvec<T, N> &rhs) {

    lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

    for (size_t i = 0; i < N - 1; i++) {

        lhs << rhs.get(i) << ", ";
    }

    return lhs << rhs.get(N - 1) << ")";
}

#define __m_impl__
#endif
