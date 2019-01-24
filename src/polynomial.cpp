
#include <m/m.h>

#include <m/polynomial.h>

m::ComplexSolutions::ComplexSolutions() noexcept
    :variableName('z'), inf(false) {}

m::ComplexSolutions::ComplexSolutions(std::unordered_set<m::comp> finiteSet) noexcept
    :variableName('z'), solutionSet(finiteSet), inf(false) {}

char m::ComplexSolutions::getVariableName() const {

    return variableName;
}

m::ComplexSolutions m::ComplexSolutions::setVariableName(char newName) {

    if ((newName >= 'A' && newName <= 'Z') || (newName >= 'a' && newName <= 'z')) {

        variableName = newName;

    } else {

        throw std::invalid_argument("m::exception: cannot use non-alphabet variable name");
    }

    return *this;
}

m::ComplexSolutions m::ComplexSolutions::empty() noexcept {

    return ComplexSolutions();
}

m::ComplexSolutions m::ComplexSolutions::finite(std::unordered_set<m::comp> finiteSet) noexcept {

    return ComplexSolutions(finiteSet);
}

m::ComplexSolutions m::ComplexSolutions::infinite() noexcept {

    ComplexSolutions result;

    result.inf = true;

    return result;
}

// TODO: Parametrize variable name or something idk
std::ostream &m::operator<<(std::ostream &stream, const m::ComplexSolutions &solutions) {

    if (solutions.inf) {

        stream << solutions.variableName << " in C";

    } else if (solutions.solutionSet.empty()) {

        stream << "no such " << solutions.variableName;

    } else {

        size_t n = solutions.solutionSet.size();
        size_t i = 0;

        for (auto root : solutions.solutionSet) {

            if (i++ == n - 1) {

                stream << solutions.variableName << " = " << root;

            } else {

                stream << solutions.variableName << " = " << root << ", or ";
            }
        }
    }

    return stream;
}

bool m::operator==(const m::PolynomialDegree &lhs, const m::PolynomialDegree &rhs) {

    if (lhs.isInfinite()) return rhs.isInfinite();

    return lhs.getValue() == rhs.getValue();
}

bool m::operator!=(const m::PolynomialDegree &lhs, const m::PolynomialDegree &rhs) {

    return !(lhs == rhs);
}

std::ostream &m::operator<<(std::ostream &lhs, const m::PolynomialDegree &rhs) {

    if (rhs.inf) return lhs << "infinity";

    return lhs << rhs.value;
}

m::PolynomialDegree::PolynomialDegree(size_t value)
    :value(value), inf(false) {}

bool m::PolynomialDegree::isInfinite() const {

    return inf;
}

size_t m::PolynomialDegree::getValue() const {

    return value;
}

m::PolynomialDegree m::PolynomialDegree::infinite() {

    PolynomialDegree result(0);
    result.inf = true;
    return result;
}

m::Polynomial::Polynomial() noexcept
    :variableName('z'), rootsValid(true), roots(ComplexSolutions::infinite()), degreeValid(true), degree(PolynomialDegree::infinite()) {}

m::Polynomial m::Polynomial::fromCoeffs(std::vector<m::comp> coeffs) noexcept {

    m::Polynomial result;

    result.coeffs = coeffs;
    result.rootsValid = false;
    result.roots = ComplexSolutions::empty();
    result.degreeValid = false;
    result.degree = 0;

    return result;
}

m::Polynomial::operator std::function<m::comp(m::comp)>() const {

    return [&] (comp z) { return this->value(z); };
}

m::comp m::Polynomial::operator()(const m::comp &z) const {

    return value(z);
}

void m::Polynomial::updateValues() {

    while (util::isZero(coeffs.back())) {

        coeffs.pop_back();
    }
}

void m::Polynomial::updateDegree() noexcept {

    degreeValid = true;

    auto l = coeffs.size();

    while (util::isZero(coeffs[l - 1])) {

        if (l == 1) {

            degree = PolynomialDegree::infinite();
            return;
        }

        l--;
    }

    degree = PolynomialDegree(l - 1);
}

std::vector<m::comp> m::Polynomial::getCoeffs() {

    updateValues();

    return coeffs;
}

std::vector<m::comp> m::Polynomial::getCoeffs() const {

    return coeffs;
}

m::Polynomial m::Polynomial::setVariableName(char newName) {

    if ((newName >= 'A' && newName <= 'Z') || (newName >= 'a' && newName <= 'z')) {

        variableName = newName;

    } else {

        throw std::invalid_argument("m::exception: cannot use non-alphabet variable name");
    }

    return *this;
}

char m::Polynomial::getVariableName() const {

    return variableName;
}

m::PolynomialDegree m::Polynomial::getDegree() {

    if (!degreeValid) updateDegree();

    return degree;
}

m::PolynomialDegree m::Polynomial::getDegree() const {

    Polynomial c = *this;

    return c.getDegree();
}

m::comp m::Polynomial::value(m::comp x) const {

    m::comp result = 0;
    m::comp v = 1;

    for (auto coeff : coeffs) {

        result += v * coeff;
        v *= x;
    }

    return result;
}

m::ComplexSolutions m::Polynomial::solve() {

    if (rootsValid) return roots;

    rootsValid = true;

    if (!degreeValid) updateDegree();

    if (degree.isInfinite()) return roots = ComplexSolutions::infinite().setVariableName(variableName);

    switch (degree.getValue()) {

        case 0: {

            // Non-zero because of is for infinite degree above
            return roots = ComplexSolutions::empty().setVariableName(variableName);
        }

        case 1: {

            return roots = ComplexSolutions::finite(-coeffs[0] / coeffs[1]).setVariableName(variableName);
        }

        case 2: {

                auto descriminant = coeffs[1] * coeffs[1] - 4.0 * coeffs[2] * coeffs[0];
                auto offset = std::sqrt(descriminant);

                auto lesser = -coeffs[1] - offset;
                auto greater = -coeffs[1] + offset;

                auto denom = 2.0 * coeffs[2];

                if (util::isEqual(lesser, greater)) return roots = ComplexSolutions::finite(lesser / denom).setVariableName(variableName);

                return roots = ComplexSolutions::finite(lesser / denom, greater / denom).setVariableName(variableName);
        }

        default: {

            // TODO: Numerical solution for higher degrees
            // TODO: Cubics and quartics
            
            return roots = ComplexSolutions::empty().setVariableName(variableName);
        }
    }
}

m::comp m::Polynomial::getCoeff(size_t index) const {

    if (index >= coeffs.size()) return comp::fromCartesian(0, 0);

    return coeffs[index];
}

void m::Polynomial::setCoeff(size_t index, const m::comp &value) {

    if (index >= coeffs.size()) {

        while (coeffs.size() < index) {

            coeffs.push_back(comp::fromCartesian(0, 0));
        }

        coeffs.push_back(value);

    } else {

        coeffs[index] = value;
    }

    rootsValid = false;
    degreeValid = false;
}

m::Polynomial m::Polynomial::interpolate(const std::vector<cvec2> &points) {

    return interpolate(points, 0, points.size() - 1);
}

m::Polynomial m::Polynomial::interpolate(const std::vector<cvec2> &points, size_t first, size_t last) {

    if (first == last) return fromCoeffs(points[first].y());

    auto leftP = interpolate(points, first, last - 1);
    auto rightP = interpolate(points, first + 1, last);

    auto leftX = points[first].x();
    auto rightX = points[last].x();

    auto leftLinear = fromCoeffs(leftX, comp(-1));
    auto rightLinear = fromCoeffs(-rightX, comp(1));

    return (rightLinear * leftP + leftLinear * rightP) / (leftX - rightX);
}

m::Polynomial &m::Polynomial::operator+=(const m::Polynomial &rhs) {

    return *this = *this + rhs;
}

m::Polynomial &m::Polynomial::operator-=(const m::Polynomial &rhs) {

    return *this = *this - rhs;
}

m::Polynomial &m::Polynomial::operator*=(const m::Polynomial &rhs) {

    return *this = *this * rhs;
}

m::Polynomial &m::Polynomial::operator+=(const m::comp &rhs) {

    return *this = *this + rhs;
}

m::Polynomial &m::Polynomial::operator-=(const m::comp &rhs) {

    return *this = *this - rhs;
}

m::Polynomial &m::Polynomial::operator*=(const m::comp &rhs) {

    return *this = *this * rhs;
}

m::Polynomial &m::Polynomial::operator/=(const m::comp &rhs) {

    return *this = *this / rhs;
}

m::Polynomial m::operator+(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    result.setCoeff(0, result.getCoeff(0) + rhs);

    return result;
}

m::Polynomial m::operator+(const m::comp &lhs, const m::Polynomial &rhs) {

    return rhs + lhs;
}

m::Polynomial m::operator+(const m::Polynomial &lhs, const m::Polynomial &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < rhs.coeffs.size(); i++) {

        result.setCoeff(i, lhs.getCoeff(i) + rhs.getCoeff(i));
    }

    return result;
}

m::Polynomial m::operator-(const m::Polynomial &rhs) {

    auto result = rhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] = -result.coeffs[i];
    }

    return result;
}

m::Polynomial m::operator-(const m::Polynomial &lhs, const m::comp &rhs) {

    return lhs + (-rhs);
}

m::Polynomial m::operator-(const m::comp &lhs, const m::Polynomial &rhs) {

    return lhs + (-rhs);
}

m::Polynomial m::operator-(const m::Polynomial &lhs, const m::Polynomial &rhs) {

    return lhs + (-rhs);
}

m::Polynomial m::operator*(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] *= rhs;
    }

    return result;
}

m::Polynomial m::operator*(const m::comp &lhs, const m::Polynomial &rhs) {

    return rhs * lhs;
}

m::Polynomial m::operator*(const m::Polynomial &lhs, const m::Polynomial &rhs) {

    // TODO: Multivariable polynomials; compare variable names

    if (lhs.getDegree().isInfinite()) return rhs;
    if (rhs.getDegree().isInfinite()) return lhs;

    Polynomial result;

    auto N = lhs.getDegree().getValue();
    auto M = rhs.getDegree().getValue();

    for (size_t i = 0; i <= N + M; i++) {

        comp c = 0;

        for (size_t j = 0; j <= i; j++) {

            size_t k = i - j;

            auto a = j > N ? 0 : lhs.getCoeff(j);
            auto b = k > M ? 0 : rhs.getCoeff(k);

            c += a * b;
        }

        result.setCoeff(i, c);
    }

    return result;
}

m::Polynomial m::operator/(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] /= rhs;
    }

    return result;
}

bool m::operator==(const Polynomial &lhs, const Polynomial &rhs) {

    auto lDeg = lhs.getDegree();
    auto rDeg = rhs.getDegree();

    if (lDeg != rDeg) return false;

    if (lDeg.isInfinite()) return true;

    auto N = lDeg.getValue();

    for (size_t i = 0; i < N; i++) {

        if (!util::isEqual(lhs.coeffs[i], rhs.coeffs[i])) return false;
    }

    return true;
}

bool m::operator!=(const Polynomial &lhs, const Polynomial &rhs) {


    return !(lhs == rhs);
}

std::ostream &m::operator<<(std::ostream &lhs, const m::Polynomial &rhs) {

    if (rhs.getDegree().isInfinite()) return lhs << "0";

    auto N = rhs.getDegree().getValue();

    bool nonZero = false;

    auto cTerm = rhs.coeffs[0];
   
    if (!util::isZero(cTerm)) {
        
        lhs << cTerm;
        nonZero = true;
    }

    if (N == 0) return lhs;

    auto lTerm = rhs.coeffs[1];

    if (!util::isZero(lTerm)) {

        if (nonZero) lhs << " + ";
        if (!util::isEqual(lTerm, comp::fromCartesian(1, 0))) lhs << lTerm;
        lhs << rhs.variableName;

        nonZero = true;
    }

    for (size_t i = 2; i < N + 1; i++) {

        auto term = rhs.coeffs[i];

        if (util::isZero(term)) continue;

        if (nonZero) lhs << " + ";
        if (!util::isEqual(term, comp::fromCartesian(1, 0))) lhs << term;
        lhs << rhs.variableName << "^" << i;

        nonZero = true;
    }

    return lhs;
}

m::Polynomial m::differentiate(const m::Polynomial &polynomial) {

    if (polynomial.getDegree().isInfinite()) return Polynomial();

    auto N = polynomial.getDegree().getValue();

    Polynomial result;

    for (size_t i = 1; i < N + 1; i++) {

        result.setCoeff(i - 1, m::comp::fromCartesian(i, 0) * polynomial.getCoeff(i));
    }

    return result;
}

m::Polynomial m::integrate(const m::Polynomial &polynomial) {

    if (polynomial.getDegree().isInfinite()) return Polynomial();

    auto N = polynomial.getDegree().getValue();

    Polynomial result;

    for (size_t i = 0; i < N + 1; i++) {

        result.setCoeff(i + 1, polynomial.getCoeff(i) / m::comp::fromCartesian(i + 1, 0)); 
    }

    return result;
}
