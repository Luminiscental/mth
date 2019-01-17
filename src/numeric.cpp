
#include <m/m.h>

#include <m/numeric.h>
#include <m/polynomial.h>

m::Polynomial lerpTowards(std::function<m::comp(m::comp)> xTransform, std::function<m::comp(size_t,m::comp)> yFunc) {

    std::vector<m::cvec2> result;
    m::comp lastX, lastY;

    for (size_t i = 2; i < 100; i++) {

        auto approachingZero = std::pow(m::comp(2), -m::comp(i));
        auto x = xTransform(approachingZero);
        auto y = yFunc(i, x);

        // If sequence is close enough to cause division by zero then we can probably break
        if (std::isnan(y.real()) || std::isnan(y.imag()) || std::isnan(x.real()) || std::isnan(x.imag())) {

            break;
        }

        result.push_back(m::cvec2(x, y));

        auto deltaX = std::abs(x - lastX);
        auto deltaY = std::abs(y - lastY);

        lastX = x;
        lastY = y;
    }

    // number of vertices to interpolate
    auto n = 5;

    auto lastIndex = result.size() - 1;
    auto startIndex = lastIndex > (n - 1) ? lastIndex - n : 0;

    return m::Polynomial::interpolate(result, startIndex, lastIndex);
}

std::function<m::comp(size_t)> accelerateConvergence(const std::function<m::comp(size_t)> &sequence) {

    // Shank transform
    return [sequence] (size_t n) {

        if (n == 0) return m::comp(0);

        auto next = sequence(n + 1);
        auto curr = sequence(n);
        auto prev = sequence(n - 1);

        auto step = next - curr;

        return next - step * step / (step - curr + prev);
    };
}

m::comp m::numeric::limit(const std::function<m::comp(size_t)> &sequence) {

    auto accelerated = accelerateConvergence(sequence);

    auto id = [] (comp z) { return z; };
    auto y = [&] (size_t index, comp x) { return accelerated(index); };

    auto pol = lerpTowards(id, y);

    return pol.getCoeff(0);
}

m::comp m::numeric::lowerLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input - small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto pol = lerpTowards(x, y);

    return pol.getCoeff(0);
}

m::comp m::numeric::upperLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input + small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto pol = lerpTowards(x, y);

    return pol.getCoeff(0);
}

m::comp m::numeric::limit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    return lowerLimit(function, input);
}

m::comp m::numeric::limitInfPos(const std::function<m::comp(m::comp)> &function) {

    auto inverted = [function] (m::comp z) { return function(z.inverse()); };
    return upperLimit(inverted, comp(0));
}

m::comp m::numeric::limitInfNeg(const std::function<m::comp(m::comp)> &function) {

    auto inverted = [function] (m::comp z) { return function(z.inverse()); };
    return lowerLimit(inverted, comp(0));
}

std::function<m::comp(m::comp)> m::differentiate(const std::function<m::comp(m::comp)> &function) {

    auto avgGradient = [&] (comp a, comp b) {

        auto dy = function(b) - function(a);
        auto dx = b - a;

        return dy / dx;
    };

    return [&] (comp x) {

        auto gradientApprox = [&] (m::comp dx) {

            return avgGradient(x, x + dx);
        };

        return numeric::limit(gradientApprox, comp(0));
    };
}

