
#include <mth/mth.h>

#include <mth/polynomial.h>

mth::ComplexSolutions::ComplexSolutions(std::unordered_set<mth::comp> finiteSet) noexcept
    :solutionSet(finiteSet) {}

char mth::ComplexSolutions::getVariableName() const {

    return variableName;
}

bool mth::ComplexSolutions::contains(const mth::comp &z) const {

    if (isInfinite()) {

        return true;

    } else {

        return solutionSet.count(z) != 0;
    }
}

bool mth::ComplexSolutions::isInfinite() const {

    return inf;
}

mth::ComplexSolutions mth::ComplexSolutions::setVariableName(char newName) {

    if ((newName >= 'A' && newName <= 'Z') || (newName >= 'a' && newName <= 'z')) {

        variableName = newName;

    } else {

        throw std::invalid_argument("mth::exception: cannot use non-alphabet variable name");
    }

    return *this;
}

mth::ComplexSolutions mth::ComplexSolutions::empty() noexcept {

    return ComplexSolutions();
}

mth::ComplexSolutions mth::ComplexSolutions::finite(std::unordered_set<mth::comp> finiteSet) noexcept {

    return ComplexSolutions(finiteSet);
}

mth::ComplexSolutions mth::ComplexSolutions::infinite() noexcept {

    ComplexSolutions result;

    result.inf = true;

    return result;
}

// TODO: Parametrize variable name or something idk
std::ostream &mth::operator<<(std::ostream &stream, const mth::ComplexSolutions &solutions) {

    if (solutions.inf) {

        stream << solutions.variableName << " in C";

    } else if (solutions.solutionSet.empty()) {

        stream << "no such " << solutions.variableName;

    } else {

        auto n = solutions.solutionSet.size();
        auto i = size_t{0};

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

bool mth::operator==(const mth::PolynomialDegree &lhs, const mth::PolynomialDegree &rhs) {

    if (lhs.isInfinite()) return rhs.isInfinite();

    return lhs.getValue() == rhs.getValue();
}

bool mth::operator!=(const mth::PolynomialDegree &lhs, const mth::PolynomialDegree &rhs) {

    return !(lhs == rhs);
}

std::ostream &mth::operator<<(std::ostream &lhs, const mth::PolynomialDegree &rhs) {

    if (rhs.inf) return lhs << "infinity";

    return lhs << rhs.value;
}

mth::PolynomialDegree::PolynomialDegree(size_t value)
    :value(value) {}

bool mth::PolynomialDegree::isInfinite() const {

    return inf;
}

size_t mth::PolynomialDegree::getValue() const {

    return value;
}

mth::PolynomialDegree mth::PolynomialDegree::infinite() {

    PolynomialDegree result(0);
    result.inf = true;
    return result;
}

mth::Polynomial mth::Polynomial::fromCoeffs(std::vector<mth::comp> coeffs) noexcept {

    mth::Polynomial result;

    result.coeffs = coeffs;

    // Lazily evaluate these; i.e. don't now
    result.rootsValid = false;
    result.roots = ComplexSolutions::empty();
    result.degreeValid = false;
    result.degree = 0;

    return result;
}

mth::Polynomial::operator std::function<mth::comp(mth::comp)>() const {

    return [&] (comp z) { return this->value(z); };
}

mth::comp mth::Polynomial::operator()(const mth::comp &z) const {

    return value(z);
}

void mth::Polynomial::updateValues() {

    while (util::isZero(coeffs.back())) {

        coeffs.pop_back();
    }
}

void mth::Polynomial::updateDegree() noexcept {

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

const std::vector<mth::comp> &mth::Polynomial::getCoeffs() {

    updateValues();

    return coeffs;
}

const std::vector<mth::comp> &mth::Polynomial::getCoeffs() const {

    return coeffs;
}

mth::Polynomial mth::Polynomial::setVariableName(char newName) {

    if ((newName >= 'A' && newName <= 'Z') || (newName >= 'a' && newName <= 'z')) {

        variableName = newName;

    } else {

        throw std::invalid_argument("mth::exception: cannot use non-alphabet variable name");
    }

    return *this;
}

char mth::Polynomial::getVariableName() const {

    return variableName;
}

mth::PolynomialDegree mth::Polynomial::getDegree() {

    if (!degreeValid) updateDegree();

    return degree;
}

mth::PolynomialDegree mth::Polynomial::getDegree() const {

    auto c = *this;

    return c.getDegree();
}

mth::comp mth::Polynomial::value(mth::comp z) const {

    auto result = mth::comp{0};
    auto v = mth::comp{1};

    for (auto coeff : coeffs) {

        result += v * coeff;
        v *= z;
    }

    return result;
}

mth::ComplexSolutions mth::Polynomial::solve() {

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

            using std::sqrt;

            auto descriminant = coeffs[1] * coeffs[1] - 4.0 * coeffs[2] * coeffs[0];
            auto offset = sqrt(descriminant);

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

mth::comp mth::Polynomial::getCoeff(size_t index) const {

    if (index >= coeffs.size()) return comp::fromCartesian(0, 0);

    return coeffs[index];
}

void mth::Polynomial::setCoeff(size_t index, const mth::comp &value) {

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

mth::Polynomial mth::Polynomial::interpolate(const std::vector<cvec2> &points) {

    return interpolate(points, 0, points.size() - 1);
}

mth::Polynomial mth::Polynomial::interpolate(const std::vector<cvec2> &points, size_t first, size_t last) {

    if (first == last) return fromCoeffs(points[first].y());

    auto leftP = interpolate(points, first, last - 1);
    auto rightP = interpolate(points, first + 1, last);

    auto leftX = points[first].x();
    auto rightX = points[last].x();

    auto leftLinear = fromCoeffs(leftX, comp(-1));
    auto rightLinear = fromCoeffs(-rightX, comp(1));

    return (rightLinear * leftP + leftLinear * rightP) / (leftX - rightX);
}

mth::Polynomial &mth::Polynomial::operator+=(const mth::Polynomial &rhs) {

    return *this = *this + rhs;
}

mth::Polynomial &mth::Polynomial::operator-=(const mth::Polynomial &rhs) {

    return *this = *this - rhs;
}

mth::Polynomial &mth::Polynomial::operator*=(const mth::Polynomial &rhs) {

    return *this = *this * rhs;
}

mth::Polynomial &mth::Polynomial::operator+=(const mth::comp &rhs) {

    return *this = *this + rhs;
}

mth::Polynomial &mth::Polynomial::operator-=(const mth::comp &rhs) {

    return *this = *this - rhs;
}

mth::Polynomial &mth::Polynomial::operator*=(const mth::comp &rhs) {

    return *this = *this * rhs;
}

mth::Polynomial &mth::Polynomial::operator/=(const mth::comp &rhs) {

    return *this = *this / rhs;
}

mth::Polynomial mth::operator+(const mth::Polynomial &lhs, const mth::comp &rhs) {

    auto result = lhs;

    result.setCoeff(0, result.getCoeff(0) + rhs);

    return result;
}

mth::Polynomial mth::operator+(const mth::comp &lhs, const mth::Polynomial &rhs) {

    return rhs + lhs;
}

mth::Polynomial mth::operator+(const mth::Polynomial &lhs, const mth::Polynomial &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < rhs.coeffs.size(); i++) {

        result.setCoeff(i, lhs.getCoeff(i) + rhs.getCoeff(i));
    }

    return result;
}

mth::Polynomial mth::operator-(const mth::Polynomial &rhs) {

    auto result = rhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] = -result.coeffs[i];
    }

    return result;
}

mth::Polynomial mth::operator-(const mth::Polynomial &lhs, const mth::comp &rhs) {

    return lhs + (-rhs);
}

mth::Polynomial mth::operator-(const mth::comp &lhs, const mth::Polynomial &rhs) {

    return lhs + (-rhs);
}

mth::Polynomial mth::operator-(const mth::Polynomial &lhs, const mth::Polynomial &rhs) {

    return lhs + (-rhs);
}

mth::Polynomial mth::operator*(const mth::Polynomial &lhs, const mth::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] *= rhs;
    }

    return result;
}

mth::Polynomial mth::operator*(const mth::comp &lhs, const mth::Polynomial &rhs) {

    return rhs * lhs;
}

mth::Polynomial mth::operator*(const mth::Polynomial &lhs, const mth::Polynomial &rhs) {

    // TODO: Multivariable polynomials; compare variable names

    if (lhs.getDegree().isInfinite()) return rhs;
    if (rhs.getDegree().isInfinite()) return lhs;

    Polynomial result;

    auto N = lhs.getDegree().getValue();
    auto M = rhs.getDegree().getValue();

    for (size_t i = 0; i <= N + M; i++) {

        auto c = comp{0};

        for (size_t j = 0; j <= i; j++) {

            auto k = size_t{i - j};

            auto a = j > N ? 0 : lhs.getCoeff(j);
            auto b = k > M ? 0 : rhs.getCoeff(k);

            c += a * b;
        }

        result.setCoeff(i, c);
    }

    return result;
}

mth::Polynomial mth::operator/(const mth::Polynomial &lhs, const mth::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] /= rhs;
    }

    return result;
}

bool mth::operator==(const Polynomial &lhs, const Polynomial &rhs) {

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

bool mth::operator!=(const Polynomial &lhs, const Polynomial &rhs) {


    return !(lhs == rhs);
}

std::ostream &mth::operator<<(std::ostream &lhs, const mth::Polynomial &rhs) {

    if (rhs.getDegree().isInfinite()) return lhs << "0";

    auto N = rhs.getDegree().getValue();

    auto nonZero = false;

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

mth::Polynomial mth::differentiate(const mth::Polynomial &polynomial) {

    if (polynomial.getDegree().isInfinite()) return Polynomial();

    auto N = polynomial.getDegree().getValue();

    Polynomial result;

    for (size_t i = 1; i < N + 1; i++) {

        result.setCoeff(i - 1, mth::comp::fromCartesian(i, 0) * polynomial.getCoeff(i));
    }

    return result;
}

mth::Polynomial mth::integrate(const mth::Polynomial &polynomial) {

    if (polynomial.getDegree().isInfinite()) return Polynomial();

    auto N = polynomial.getDegree().getValue();

    Polynomial result;

    for (size_t i = 0; i < N + 1; i++) {

        result.setCoeff(i + 1, polynomial.getCoeff(i) / mth::comp::fromCartesian(i + 1, 0)); 
    }

    return result;
}

