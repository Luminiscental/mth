#ifdef __m_impl__
#undef __m_impl__

size_t std::hash<m::comp>::operator()(const m::comp &z) const {

    auto r = z.real();
    auto i = z.imag();

    return hash<decltype(r)>()(r) ^ hash<decltype(i)>()(i);
}

inline m::ComplexSolutions::ComplexSolutions() noexcept
    :inf(false) {}

inline m::ComplexSolutions::ComplexSolutions(std::unordered_set<m::comp> finiteSet) noexcept
    :solutionSet(finiteSet), inf(false) {}

inline m::ComplexSolutions m::ComplexSolutions::empty() noexcept {

    return ComplexSolutions();
}

inline m::ComplexSolutions m::ComplexSolutions::finite(std::unordered_set<m::comp> finiteSet) noexcept {

    return ComplexSolutions(finiteSet);
}

template <typename ...Q>
m::ComplexSolutions m::ComplexSolutions::finite(Q... args) noexcept {

    return ComplexSolutions(std::unordered_set<m::comp>({args...}));
}

inline m::ComplexSolutions m::ComplexSolutions::infinite() noexcept {

    ComplexSolutions result;

    result.inf = true;

    return result;
}

// TODO: Parametrize variable name or something idk
auto &m::operator<<(std::ostream &stream, const m::ComplexSolutions &solutions) {

    stream << std::fixed << std::setprecision(m_PRECISION);

    if (solutions.inf) {

        stream << "z in C";

    } else if (solutions.solutionSet.empty()) {

        stream << "no such z";

    } else {

        size_t n = solutions.solutionSet.size();
        size_t i = 0;

        for (auto root : solutions.solutionSet) {

            if (i++ == n - 1) {

                stream << "z = " << root;

            } else {

                stream << "z = " << root << ", or ";
            }
        }
    }

    return stream;
}

auto m::operator==(const m::PolynomialDegree &lhs, const m::PolynomialDegree &rhs) {

    if (lhs.isInfinite()) return rhs.isInfinite();

    return lhs.getValue() == rhs.getValue();
}

auto m::operator!=(const m::PolynomialDegree &lhs, const m::PolynomialDegree &rhs) {

    return !(lhs == rhs);
}

auto &m::operator<<(std::ostream &lhs, const m::PolynomialDegree &rhs) {

    if (rhs.inf) return lhs << "infinity";

    return lhs << rhs.value;
}

inline m::PolynomialDegree::PolynomialDegree(size_t value)
    :value(value), inf(false) {}

inline bool m::PolynomialDegree::isInfinite() const {

    return inf;
}

inline size_t m::PolynomialDegree::getValue() const {

    return value;
}

inline m::PolynomialDegree m::PolynomialDegree::infinite() {

    PolynomialDegree result(0);
    result.inf = true;
    return result;
}

m::Polynomial::Polynomial() noexcept
    :rootsValid(true), roots(ComplexSolutions::infinite()), degreeValid(true), degree(PolynomialDegree::infinite()) {}

m::Polynomial::Polynomial(std::vector<m::comp> coeffs) noexcept
    :coeffs(coeffs), rootsValid(false), roots(ComplexSolutions::empty()), degreeValid(false), degree(0) {}

template <typename ...Q>
m::Polynomial::Polynomial(Q... args) noexcept
    :Polynomial(std::vector<comp> {static_cast<comp>(args)...}) {}

m::Polynomial::operator std::function<comp(comp)>() const {

    return [&] (comp z) { return this->value(z); };
}

void m::Polynomial::updateValues() {

    while (util::checkZero(coeffs.back())) {

        coeffs.pop_back();
    }
}

void m::Polynomial::updateDegree() noexcept {

    degreeValid = true;

    auto l = coeffs.size();

    while (util::checkZero(coeffs[l - 1])) {

        if (l == 1) {

            degree = PolynomialDegree::infinite();
            return;
        }

        l--;
    }

    degree = PolynomialDegree(l - 1);
}

auto m::Polynomial::getCoeffs() {

    updateValues();

    return coeffs;
}

auto m::Polynomial::getCoeffs() const {

    return coeffs;
}

auto m::Polynomial::getDegree() {

    if (!degreeValid) updateDegree();

    return degree;
}

auto m::Polynomial::getDegree() const {

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

auto m::Polynomial::solve() {

    if (rootsValid) return roots;

    rootsValid = true;

    if (!degreeValid) updateDegree();

    if (degree.isInfinite()) return ComplexSolutions::infinite();

    switch (degree.getValue()) {

        case 0: {

            // NOTE: non-zero because of check for infinite degree above
            return ComplexSolutions::empty();
        }

        case 1: {

            return roots = ComplexSolutions::finite(-coeffs[0] / coeffs[1]);
        }

        case 2: {

                auto descriminant = coeffs[1] * coeffs[1] - 4.0f * coeffs[2] * coeffs[0];
                auto offset = std::sqrt(descriminant);

                auto lesser = -coeffs[1] - offset;
                auto greater = -coeffs[1] + offset;

                auto denom = 2.0f * coeffs[2];

                if (util::checkEqual(lesser, greater)) return roots = ComplexSolutions::finite(lesser / denom);

                return roots = ComplexSolutions::finite(lesser / denom, greater / denom);
        }

        default: {

            // TODO: Numerical solution for higher degrees
            // TODO: Cubics and quartics
            
            return roots = ComplexSolutions::empty();
        }
    }
}

auto m::Polynomial::getCoeff(size_t index) const {

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

auto m::operator+(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    result.setCoeff(0, result.getCoeff(0) + rhs);

    return result;
}

auto m::operator+(const m::comp &lhs, const m::Polynomial &rhs) {

    return rhs + lhs;
}

auto m::operator+(const m::Polynomial &lhs, const m::Polynomial &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < rhs.coeffs.size(); i++) {

        result.setCoeff(i, lhs.getCoeff(i) + rhs.getCoeff(i));
    }
}

auto m::operator-(const m::Polynomial &rhs) {

    auto result = rhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] = -result.coeffs[i];
    }

    return result;
}

auto m::operator-(const m::Polynomial &lhs, const m::comp &rhs) {

    return lhs + (-rhs);
}

auto m::operator-(const m::comp &lhs, const m::Polynomial &rhs) {

    return lhs + (-rhs);
}

auto m::operator-(const m::Polynomial &lhs, const m::Polynomial &rhs) {

    return lhs + (-rhs);
}

auto m::operator*(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] *= rhs;
    }

    return result;
}

auto m::operator*(const m::comp &lhs, const m::Polynomial &rhs) {

    return rhs * lhs;
}

auto m::operator*(const m::Polynomial &lhs, const m::Polynomial &rhs) {

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

auto m::operator/(const m::Polynomial &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < result.coeffs.size(); i++) {

        result.coeffs[i] /= rhs;
    }

    return result;
}

auto m::operator==(const Polynomial &lhs, const Polynomial &rhs) {

    auto lDeg = lhs.getDegree();
    auto rDeg = rhs.getDegree();

    if (lDeg != rDeg) return false;

    if (lDeg.isInfinite()) return true;

    auto N = lDeg.getValue();

    for (size_t i = 0; i < N; i++) {

        if (!util::checkEqual(lhs.coeffs[i], rhs.coeffs[i])) return false;
    }

    return true;
}

auto m::operator!=(const Polynomial &lhs, const Polynomial &rhs) {


    return !(lhs == rhs);
}

auto &m::operator<<(std::ostream &lhs, const m::Polynomial &rhs) {

    lhs << std::fixed << std::setprecision(m_PRECISION);

    if (rhs.getDegree().isInfinite()) return lhs << "0";

    auto N = rhs.getDegree().getValue();

    bool nonZero = false;

    auto cTerm = rhs.coeffs[0];
   
    if (!util::checkZero(cTerm)) {
        
        lhs << cTerm;
        nonZero = true;
    }

    if (N == 0) return lhs;

    auto lTerm = rhs.coeffs[1];

    if (!util::checkZero(lTerm)) {

        if (nonZero) lhs << " + ";
        if (!util::checkEqual(lTerm, comp::fromCartesian(1, 0))) lhs << lTerm;
        lhs << "z";

        nonZero = true;
    }

    for (size_t i = 2; i < N + 1; i++) {

        auto term = rhs.coeffs[i];

        if (util::checkZero(term)) continue;

        if (nonZero) lhs << " + ";
        if (!util::checkEqual(term, comp::fromCartesian(1, 0))) lhs << term;
        lhs << "z^" << i;

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

#define __m_impl__
#endif
