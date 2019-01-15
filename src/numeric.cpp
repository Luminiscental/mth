
#include <iostream>

#include <m/m.h>

#include <m/numeric.h>
#include <m/polynomial.h>

std::vector<m::cvec2> sampleFunction(std::function<m::comp(m::comp)> xTransform, std::function<m::comp(size_t,m::comp)> yFunc) {

    std::vector<m::cvec2> result;
    m::comp last;

    // TODO: Parametrize these numbers / choose intelligently
    for (size_t i = 1; i < 999; i++) {

        auto x = xTransform(std::pow(m::comp(2), -m::comp(i)));
        if (std::abs(x - last) < 0.00001) break;

        auto y = yFunc(i, x);
        result.push_back(m::cvec2(x, y));

        last = x;
    }

    return result;
}

m::comp m::numeric::limit(const std::function<m::comp(size_t)> &sequence) {

    auto id = [] (comp z) { return z; };
    auto y = [&] (size_t index, comp x) { return sequence(index); };

    auto testValues = sampleFunction(id, y);

    return Polynomial::interpolate(testValues).getCoeff(0);
}

m::comp m::numeric::lowerLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input - small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto testValues = sampleFunction(x, y);

    return Polynomial::interpolate(testValues).getCoeff(0);
}

m::comp m::numeric::upperLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input + small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto testValues = sampleFunction(x, y);

    return Polynomial::interpolate(testValues).getCoeff(0);
}

m::comp m::numeric::limit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    return lowerLimit(function, input);
}

m::comp m::numeric::limitInfPos(const std::function<m::comp(m::comp)> &function) {

    auto sequenceEquiv = [&] (size_t n) {

        return function(comp(n));
    };

    return limit(sequenceEquiv);
}

m::comp m::numeric::limitInfNeg(const std::function<m::comp(m::comp)> &function) {

    auto sequenceEquiv = [&] (size_t n) {

        return function(-comp(n));
    };

    return limit(sequenceEquiv);
}

std::function<m::comp(m::comp)> m::numeric::differentiate(const std::function<m::comp(m::comp)> &function) {

    auto avgGradient = [&] (comp a, comp b) {

        auto dy = function(b) - function(a);
        auto dx = b - a;

        return dy / dx;
    };

    return [&] (comp x) {

        auto gradientApprox = [&] (m::comp dx) {

            return avgGradient(x, x + dx);
        };

        return limit(gradientApprox, comp(0));
    };
}

