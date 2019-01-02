#ifdef __m_impl__
#undef __m_impl__

#include <cmath>
#include <iomanip>

#define BINDING(name, value)    template <typename T> const auto &m::tquat<T>:: name () const { return value ; } \
                                template <typename T> auto &m::tquat<T>:: name () { return value ; }

        BINDING(real, w)
        BINDING(imaginary, ijk)
        BINDING(i, ijk.x())
        BINDING(j, ijk.y())
        BINDING(k, ijk.z())

#undef BINDING

template <typename T>
auto m::tquat<T>::magnSqr() const noexcept {

    return w * w + ijk.magnSqr();
}

template <typename T>
auto m::tquat<T>::magn() const noexcept {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
auto m::tquat<T>::inverse() const {

    auto ls = magnSqr();

    if (util::checkZero(ls)) throw std::invalid_argument("zero quaternion has no inverse");

    return conjugate() / ls;
}

template <typename T>
auto m::tquat<T>::unit() const {

    auto l = static_cast<T>(magn());

    if (util::checkZero(l)) throw std::invalid_argument("zero quaternion has no unit equivalent");

    return *this / l;
}

template <typename T>
auto m::tquat<T>::rotate(const m::tvec<T, 3> &vector) const {

    return (*this * tquat<T>(vector) * inverse()).imaginary();
}

template <typename T>
auto m::tquat<T>::identity() {
    
    return tquat(); 
}

template <typename T>
auto m::tquat<T>::rotation(T angle, const m::tvec<T, 3> &axis) {

    auto halfAngle = static_cast<double>(angle) / 2.0;

    auto c = static_cast<T>(std::cos(halfAngle));
    auto s = static_cast<T>(std::sin(halfAngle));

    return c + s * tquat<T>(axis);
}

template <typename T>
auto &m::tquat<T>::operator+=(const m::tquat<T> &rhs) {

    w += rhs.w;
    ijk += rhs.ijk;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator+=(T rhs) {

    w += rhs;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator-=(const m::tquat<T> &rhs) {

    w -= rhs.w;
    ijk -= rhs.ijk;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator-=(T rhs) {

    w -= rhs;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator*=(const m::tquat<T> &rhs) {

    *this = *this * rhs;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator*=(T rhs) {

    w *= rhs;
    ijk *= rhs;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator/=(const m::tquat<T> &rhs) {

    *this = *this / rhs;

    return *this;
}

template <typename T>
auto &m::tquat<T>::operator/=(T rhs) {

    w /= rhs;
    ijk /= rhs;

    return *this;
}

template <typename T>
auto m::operator+(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result += rhs;
}

template <typename T>
auto m::operator+(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result += rhs;
}

template <typename T>
auto m::operator+(T lhs, const m::tquat<T> &rhs) {

    return rhs + lhs;
}

template <typename T>
auto m::operator-(const m::tquat<T> &rhs) {

    return tquat<T>(-rhs.w, -rhs.ijk);
}

template <typename T>
auto m::operator-(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator-(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator-(T lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
auto m::operator*(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    return tquat<T>(lhs.real() * rhs.i() + lhs.i() * rhs.real() + lhs.j() * rhs.k() - lhs.k() * rhs.j(),
                    lhs.real() * rhs.j() - lhs.i() * rhs.k() + lhs.j() * rhs.real() + lhs.k() * rhs.i(),
                    lhs.real() * rhs.k() + lhs.i() * rhs.j() - lhs.j() * rhs.i() + lhs.k() * rhs.real(),
                    lhs.real() * rhs.real() - lhs.i() * rhs.i() - lhs.j() * rhs.j() - lhs.k() * rhs.k());
}

template <typename T>
auto m::operator*(T lhs, const m::tquat<T> &rhs) {

    tquat<T> result = rhs;

    return result *= lhs;
}

template <typename T>
auto m::operator*(const m::tquat<T> &lhs, T rhs) {

    return rhs * lhs;
}

template <typename T>
auto m::operator/(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
auto m::operator/(T lhs, const m::tquat<T> &rhs) {

    return tquat<T>(lhs) / rhs;
}

template <typename T>
auto m::operator/(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result /= rhs;
}

template <typename T>
auto &m::operator<<(std::ostream &lhs, const m::tquat<T> &rhs) {

    bool nonZero = false;

    lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

    if (!util::checkZero(rhs.real())) {

        lhs << rhs.real();
        nonZero = true;
    }

    if (!util::checkZero(rhs.i())) {

        if (nonZero) {

            lhs << (rhs.i() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.i());

        } else {

            lhs << rhs.i();
            nonZero = true;
        }

        lhs << "i";
    }

    if (!util::checkZero(rhs.j())) {
        
        if (nonZero) {

            lhs << (rhs.j() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.j());

        } else {

            lhs << rhs.j();
            nonZero = true;
        }

        lhs << "j";
    }

    if (!util::checkZero(rhs.k())) {
        
        if (nonZero) {

            lhs << (rhs.k() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.k());

        } else {

            lhs << rhs.k();
        }

        lhs << "k";
    }

    if (!nonZero) lhs << "0";

    return lhs << ")";
}

#define __m_impl__
#endif
