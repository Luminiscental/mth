
#include <m/m.h>

#include <m/numeric.h>
#include <m/polynomial.h>

m::Polynomial lerpTowards(std::function<m::comp(m::comp)> xTransform, std::function<m::comp(size_t,m::comp)> yFunc) {

    using std::pow;
    using std::abs;

    std::vector<m::cvec2> result;
    m::comp lastX, lastY;

    for (size_t i = 2; i < 100; i++) {

        // TODO: Parametrize this sequence or choose it contextually
        auto approachingZero = pow(m::comp(2), -m::comp(i));
        auto x = xTransform(approachingZero);
        auto y = yFunc(i, x);

        // If sequence is close enough to cause division by zero then we can probably break
        if (std::isnan(y.real()) || std::isnan(y.imag()) || std::isnan(x.real()) || std::isnan(x.imag())) {

            break;
        }

        result.push_back(m::cvec2(x, y));

        auto deltaX = abs(x - lastX);
        auto deltaY = abs(y - lastY);

        lastX = x;
        lastY = y;
    }

    // number of vertices to interpolate
    auto n = 6;

    auto lastIndex = result.size() - 1;
    auto startIndex = lastIndex > (n - 1) ? lastIndex - n : 0;

    return m::Polynomial::interpolate(result, startIndex, lastIndex);
}

std::function<m::comp(size_t)> shankTransform(const std::function<m::comp(size_t)> &partialSum, const std::function<m::comp(size_t)> &sequence) {

    return [partialSum, sequence] (size_t n) {

        if (n == 0) return sequence(0);

        auto nextValue = sequence(n + 1);
        auto currentValue = sequence(n);
        auto nextSum = nextValue + currentValue + partialSum(n - 1);

        return nextSum - nextValue * nextValue / (nextValue - currentValue);
    };
}

std::function<m::comp(size_t)> aitkenTransform(const std::function<m::comp(size_t)> &sequence) {

    return [sequence] (size_t n) {

        if (n == 0) return m::comp(0);

        auto next = sequence(n + 1);
        auto curr = sequence(n);
        auto prev = sequence(n - 1);

        auto step = next - curr;

        return next - step * step / (step - curr + prev);
    };
}

m::comp m::limit(const std::function<m::comp(size_t)> &sequence) {

    auto accelerated = aitkenTransform(sequence);

    auto id = [] (comp z) { return z; };
    auto y = [&] (size_t index, comp x) { return accelerated(index); };

    auto pol = lerpTowards(id, y);

    return pol.getCoeff(0);
}

m::comp m::seriesLimit(const std::function<m::comp(size_t)> &partialSum, const std::function<m::comp(size_t)> &sequence) {

    auto accelerated = shankTransform(partialSum, sequence);

    auto id = [] (comp z) { return z; };
    auto y = [&] (size_t index, comp x) { return accelerated(index); };

    auto pol = lerpTowards(id, y);

    return pol.getCoeff(0);
}

m::comp m::lowerLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input - small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto pol = lerpTowards(x, y);

    return pol.getCoeff(0);
}

m::comp m::upperLimit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    auto x = [&] (comp small) { return input + small; };
    auto y = [&] (size_t index, comp x) { return function(x); };

    auto pol = lerpTowards(x, y);

    return pol.getCoeff(0);
}

m::comp m::limit(const std::function<m::comp(m::comp)> &function, m::comp input) {

    return lowerLimit(function, input);
}

m::comp m::limitInfPos(const std::function<m::comp(m::comp)> &function) {

    auto inverted = [function] (m::comp z) { return function(z.inverse()); };
    return upperLimit(inverted, comp(0));
}

m::comp m::limitInfNeg(const std::function<m::comp(m::comp)> &function) {

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

        return limit(gradientApprox, comp(0));
    };
}

