
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

#define BINDING(name, index)    template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > index), int>::type> \
                                const T &m::tvec<T, N>:: name () const { return this->get(index); } \
                                template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > index), int>::type> \
                                      T &m::tvec<T, N>:: name ()       { return this->get(index); }

        BINDING(x, 0)
        BINDING(y, 1)
        BINDING(z, 2)
        BINDING(w, 3)

#undef BINDING

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 1), int>::type>
m::tvec<T, 2> m::tvec<T, N>::xy() const {

    return tvec<T, 2>(x(), y());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 2), int>::type>
m::tvec<T, 2> m::tvec<T, N>::yz() const {

    return tvec<T, 2>(y(), z());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 3), int>::type>
m::tvec<T, 2>  m::tvec<T, N>::zw() const {

    return tvec<T, 2>(z(), w());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 2), int>::type>
m::tvec<T, 3>  m::tvec<T, N>::xyz() const {

    return tvec<T, 3>(x(), y(), z());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 3), int>::type>
m::tvec<T, 3>  m::tvec<T, N>::yzw() const {

    return tvec<T, 3>(y(), z(), w());
}

template <typename T, size_t N> template <size_t n, typename std::enable_if<(n == N) && (n > 3), int>::type>
m::tvec<T, 4> m::tvec<T, N>::xyzw() const {

    return tvec<T, 4>(x(), y(), z(), w());
}

template <typename T, size_t N>
T m::tvec<T, N>::magnSqr() const noexcept {

    auto result = static_cast<T>(0);

    for (const auto &value : *this) {

        result += value * value;
    }

    return result;
}

template <typename T, size_t N>
double m::tvec<T, N>::magn() const noexcept {

    auto ls = static_cast<double>(magnSqr());

    if (util::isZero(ls)) return 0.0;

    return std::sqrt(ls);
}

template <typename T, size_t N>
T m::tvec<T, N>::dot(const m::tvec<T, N> &rhs) const {

    T result = 0;

    for (size_t i = 0; i < N; i++) {

        result += this->get(i) * rhs.get(i);
    }

    return result;
}

template <typename T, size_t N>
m::tmat<T, N, N> m::tvec<T, N>::outerProduct(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) {

    auto lMat = m::tmat<T, 1, N>(lhs);
    auto rMatT = m::tmat<T, N, 1>(rhs);

    return lMat * rMatT;
}

template <typename T, size_t N>
m::tvec<T, N> m::tvec<T, N>::unit() const {

    auto l = static_cast<T>(magn());

    if (util::isZero(l)) throw std::invalid_argument("m::exception: unit() called on zero vector");

    return *this / l;
}

template <typename T, size_t N>
T m::tvec<T, N>::dot(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) {

    return lhs.dot(rhs);
}

template <typename T, size_t N>
m::tvec<T, N> &m::tvec<T, N>::operator+=(const m::tvec<T, N> &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) += rhs.get(i);
    }

    return *this;
}

template <typename T, size_t N>
m::tvec<T, N> &m::tvec<T, N>::operator-=(const m::tvec<T, N> &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) -= rhs.get(i);
    }

    return *this;
}

template <typename T, size_t N>
m::tvec<T, N> &m::tvec<T, N>::operator*=(const T &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) *= rhs;
    }

    return *this;
}

template <typename T, size_t N>
m::tvec<T, N> &m::tvec<T, N>::operator/=(const T &rhs) {

    for (size_t i = 0; i < N; i++) {

        this->get(i) /= rhs;
    }

    return *this;
}

template <typename T>
m::tvec<T, 3> m::vec::cross(const m::tvec<T, 3> &lhs, const m::tvec<T, 3> &rhs) {

    return tvec<T, 3>(lhs.get(1) * rhs.get(2) - lhs.get(2) * rhs.get(1),
                      lhs.get(2) * rhs.get(0) - lhs.get(0) * rhs.get(2),
                      lhs.get(0) * rhs.get(1) - lhs.get(1) * rhs.get(0));
}

template <typename T>
T m::vec::det(const m::tvec<T, 2> &lhs, const m::tvec<T, 2> &rhs) {

    return lhs.get(0) * rhs.get(1) - lhs.get(1) * rhs.get(0); // tmat<T, 2, 2>(lhs, rhs).det()
}

template <typename T, size_t N>
m::tvec<T, N> m::operator+(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) noexcept {

    auto result = lhs;

    return result += rhs;
}

template <typename T, size_t N>
m::tvec<T, N> m::operator-(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) noexcept {

    auto result = lhs;

    return result -= rhs;
}

template <typename T, size_t N>
m::tvec<T, N> m::operator-(const m::tvec<T, N> &rhs) noexcept {

    auto result = rhs;

    for (size_t i = 0; i < N; i++) {

        result.get(i) = -rhs.get(i);
    }

    return result;
}

template <typename T, size_t N>
m::tvec<T, N> m::operator*(const T &lhs, const m::tvec<T, N> &rhs) noexcept {

    auto result = rhs;

    return result *= lhs;
}

template <typename T, size_t N>
m::tvec<T, N> m::operator*(const m::tvec<T, N> &lhs, const T &rhs) noexcept {

    return rhs * lhs;
}

template <typename T, size_t N>
m::tvec<T, N> m::operator/(const m::tvec<T, N> &lhs, const T &rhs) noexcept {

    tvec<T, N> result = lhs;

    return result /= rhs;
}

template <typename T, size_t N>
bool m::operator==(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) noexcept {

    for (size_t i = 0; i < N; i++) {

        if (!util::isEqual(lhs.get(i), rhs.get(i))) return false;
    }

    return true;
}

template <typename T, size_t N>
bool m::operator!=(const m::tvec<T, N> &lhs, const m::tvec<T, N> &rhs) noexcept {

    return !(lhs == rhs);
}

template <typename T, size_t N>
std::ostream &m::operator<<(std::ostream &lhs, const m::tvec<T, N> &rhs) {

    lhs << "(";

    for (size_t i = 0; i < N - 1; i++) {

        lhs << rhs.get(i) << ", ";
    }

    return lhs << rhs.get(N - 1) << ")";
}

template <typename T, size_t N>
double std::abs(const m::tvec<T, N> &x) noexcept {

    return x.magn();
}

template <typename T, size_t N>
size_t std::hash<m::tvec<T, N>>::operator()(const m::tvec<T, N> &x) {

    size_t result = 0;

    for (size_t i = 0; i < N; i++) {

        result ^= hash<T>()(x.get(i));
    }

    return result;
}
