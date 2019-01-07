
#include <m/m.h>

#include <m/comp.h>
#include <m/vec.h>

double std::abs(const m::comp &z) {

    return z.abs();
}

m::comp std::sqrt(const m::comp &z) {

    auto p = z.asPolar();

    p.x() = sqrt(p.x());
    p.y() /= 2;

    return m::comp::fromPolar(p);
}

m::comp std::exp(const m::comp &z) {

    auto c = std::cos(z.imag());
    auto s = std::sin(z.imag());

    return exp(z.real()) * (c + m::i<float> * s);
}

m::comp std::cos(const m::comp &z) {

    auto exponent = m::i<float> * z;

    return (exp(exponent) + exp(-exponent)) / 2.0f;
}

m::comp std::sin(const m::comp &z) {

    auto exponent = m::i<float> * z;

    return (exp(exponent) - exp(-exponent)) / (2.0f * m::i<float>);
}

size_t std::hash<m::comp>::operator()(const m::comp &z) const {

    auto r = z.real();
    auto i = z.imag();

    return hash<decltype(r)>()(r) ^ hash<decltype(i)>()(i);
}
