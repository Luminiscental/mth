
#include <cmath>
#include <iomanip>

#define BINDING(name, value)    template <typename T> const auto &m::tquat<T>:: name () const { return value ; } \
                                template <typename T>       auto &m::tquat<T>:: name ()       { return value ; }

        BINDING(real, w)
        BINDING(imaginary, ijk)
        BINDING(i, ijk.x())
        BINDING(j, ijk.y())
        BINDING(k, ijk.z())

#undef BINDING

template <typename T>
T m::tquat<T>::magnSqr() const noexcept {

    return w * w + ijk.magnSqr();
}

template <typename T>
double m::tquat<T>::magn() const noexcept {

    return std::sqrt(static_cast<double>(magnSqr()));
}

template <typename T>
m::tquat<T> m::tquat<T>::inverse() const {

    auto ls = magnSqr();

    if (util::isZero(ls)) throw std::invalid_argument("zero quaternion has no inverse");

    return conjugate() / ls;
}

template <typename T>
m::tquat<T> m::tquat<T>::unit() const {

    auto l = static_cast<T>(magn());

    if (util::isZero(l)) throw std::invalid_argument("zero quaternion has no unit equivalent");

    return *this / l;
}

template <typename T>
m::tvec<T, 3> m::tquat<T>::rotate(const m::tvec<T, 3> &vector) const {

    return (*this * tquat<T>(vector) * inverse()).imaginary();
}

template <typename T>
m::tquat<T> m::tquat<T>::identity() {
    
    return tquat(); 
}

// Rotations are represented as cos(angle/2) + sin(angle/2) * axis
template <typename T>
m::tquat<T> m::tquat<T>::rotation(T angle, const m::tvec<T, 3> &axis) {

    auto halfAngle = static_cast<double>(angle) / 2.0;

    auto c = static_cast<T>(std::cos(halfAngle));
    auto s = static_cast<T>(std::sin(halfAngle));

    return c + s * tquat<T>(axis);
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator+=(const m::tquat<T> &rhs) {

    w += rhs.w;
    ijk += rhs.ijk;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator+=(T rhs) {

    w += rhs;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator-=(const m::tquat<T> &rhs) {

    w -= rhs.w;
    ijk -= rhs.ijk;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator-=(T rhs) {

    w -= rhs;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator*=(const m::tquat<T> &rhs) {

    *this = *this * rhs;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator*=(T rhs) {

    w *= rhs;
    ijk *= rhs;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator/=(const m::tquat<T> &rhs) {

    *this = *this / rhs;

    return *this;
}

template <typename T>
m::tquat<T> &m::tquat<T>::operator/=(T rhs) {

    w /= rhs;
    ijk /= rhs;

    return *this;
}

template <typename T>
m::tquat<T> m::operator+(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result += rhs;
}

template <typename T>
m::tquat<T> m::operator+(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result += rhs;
}

template <typename T>
m::tquat<T> m::operator+(T lhs, const m::tquat<T> &rhs) {

    return rhs + lhs;
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &rhs) {

    return tquat<T>(-rhs.w, -rhs.ijk);
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tquat<T> m::operator-(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tquat<T> m::operator-(T lhs, const m::tquat<T> &rhs) {

    tquat<T> result = lhs;

    return result -= rhs;
}

template <typename T>
m::tquat<T> m::operator*(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    return tquat<T>(lhs.real() * rhs.real() - lhs.i() * rhs.i() - lhs.j() * rhs.j() - lhs.k() * rhs.k(),  // real
                    lhs.real() * rhs.i() + lhs.i() * rhs.real() + lhs.j() * rhs.k() - lhs.k() * rhs.j(),  // i
                    lhs.real() * rhs.j() - lhs.i() * rhs.k() + lhs.j() * rhs.real() + lhs.k() * rhs.i(),  // j
                    lhs.real() * rhs.k() + lhs.i() * rhs.j() - lhs.j() * rhs.i() + lhs.k() * rhs.real()); // k
}

template <typename T>
m::tquat<T> m::operator*(T lhs, const m::tquat<T> &rhs) {

    tquat<T> result = rhs;

    return result *= lhs;
}

template <typename T>
m::tquat<T> m::operator*(const m::tquat<T> &lhs, T rhs) {

    return rhs * lhs;
}

template <typename T>
m::tquat<T> m::operator/(const m::tquat<T> &lhs, const m::tquat<T> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T>
m::tquat<T> m::operator/(T lhs, const m::tquat<T> &rhs) {

    return tquat<T>(lhs) / rhs;
}

template <typename T>
m::tquat<T> m::operator/(const m::tquat<T> &lhs, T rhs) {

    tquat<T> result = lhs;

    return result /= rhs;
}

template <typename T>
bool m::operator==(const tquat<T> &lhs, const tquat<T> &rhs) {

    return util::isEqual(lhs.real(), rhs.real()) && lhs.imag() == rhs.imag();
}

template <typename T>
bool m::operator==(const T &lhs, const tquat<T> &rhs) {

    return util::isEqual(lhs, rhs.real()) && rhs.imag() == tvec<T, 3>(0, 0, 0);
}

template <typename T>
bool m::operator==(const tquat<T> &lhs, const T &rhs) {

    return rhs == lhs;
}

template <typename T>
bool m::operator!=(const tquat<T> &lhs, const tquat<T> &rhs) {

    return !(lhs == rhs);
}

template <typename T>
bool m::operator!=(const T &lhs, const tquat<T> &rhs) {

    return !(lhs == rhs);
}

template <typename T>
bool m::operator!=(const tquat<T> &lhs, const T &rhs) {

    return !(lhs == rhs);
}

// TODO: Make this more consistent with m::comp (i.e. bools up-front and only use brackets when needed)
template <typename T>
std::ostream &m::operator<<(std::ostream &lhs, const m::tquat<T> &rhs) {

    bool nonZero = false;

    lhs << std::fixed << std::setprecision(m_PRECISION) << "(";

    if (!util::isZero(rhs.real())) {

        lhs << rhs.real();
        nonZero = true;
    }

    if (!util::isZero(rhs.i())) {

        if (nonZero) {

            lhs << (rhs.i() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.i());

        } else {

            lhs << rhs.i();
            nonZero = true;
        }

        lhs << "i";
    }

    if (!util::isZero(rhs.j())) {
        
        if (nonZero) {

            lhs << (rhs.j() > 0 ? " + " : " - ");
            lhs << std::abs(rhs.j());

        } else {

            lhs << rhs.j();
            nonZero = true;
        }

        lhs << "j";
    }

    if (!util::isZero(rhs.k())) {
        
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
