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

auto &m::operator<<(std::ostream &lhs, const PolynomialDegree &rhs) {

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

m::Polynomial<0>::Polynomial() noexcept
    :rootsValid(false), roots(ComplexSolutions::empty()) {

    coeff = comp::fromCartesian(0, 0);
}

template <size_t M, typename std::enable_if<(M > 0), int>::type>
m::Polynomial<0>::Polynomial(std::array<m::comp, M> coeffs) noexcept
    :rootsValid(false), roots(ComplexSolutions::empty()) {

    coeff = coeffs[0];
}

m::Polynomial<0>::Polynomial(m::comp value) noexcept
    :coeff(value), rootsValid(false), roots(ComplexSolutions::empty()) {}

auto m::Polynomial<0>::actualDegree() const {

    if (util::checkZero(coeff)) return PolynomialDegree::infinite();

    return PolynomialDegree(0);
}

auto m::Polynomial<0>::value(m::comp x) {

    return coeff;
}

auto m::Polynomial<0>::solve() {

    if (rootsValid) return roots;

    rootsValid = true;

    if (util::checkZero(coeff)) return roots = ComplexSolutions::infinite();

    return roots = ComplexSolutions::empty();
}

auto m::Polynomial<0>::getCoeff(size_t index) const {

    if (index > 0) return comp::fromCartesian(0, 0);

    return coeff;
}

void m::Polynomial<0>::setCoeff(size_t index, const comp &value) {

    if (index > 0) throw std::invalid_argument("m::exception: cannot extend polynomial");

    coeff = value;
    rootsValid = false;
}

template <size_t N>
m::Polynomial<N>::Polynomial() noexcept
    :rootsValid(false), roots(ComplexSolutions::empty()) {

    coeffs.fill(comp::fromCartesian(0, 0));
}

template <size_t N> template <size_t M, typename std::enable_if<(M > N), int>::type>
m::Polynomial<N>::Polynomial(std::array<m::comp, M> coeffs) noexcept
    :rootsValid(false), roots(ComplexSolutions::empty()) {

    for (size_t i = 0; i < N + 1; i++) this->coeffs[i] = coeffs[i];
}

template <size_t N> template <typename ...Q, typename std::enable_if<sizeof...(Q) == N + 1, int>::type>
m::Polynomial<N>::Polynomial(Q... args) noexcept
    :Polynomial(std::array<m::comp, N + 1> {args...}) {}

template <size_t N>
std::array<m::comp, N + 1> m::Polynomial<N>::getCoeffs() const {

    return coeffs;
}

template <size_t N>
auto m::Polynomial<N>::actualDegree() const {

    if (util::checkZero(coeffs[N])) return Polynomial<N - 1>(coeffs).actualDegree();

    return PolynomialDegree(N);
}

template <size_t N>
auto m::Polynomial<N>::value(m::comp x) {

    m::comp result = 0;
    m::comp v = 1;

    for (size_t i = 0; i < N + 1; i++) {

        result += v * coeffs[i];
        v *= x;
    }

    return result;
}

template <size_t N>
auto m::Polynomial<N>::solve() {

    if (rootsValid) return roots;

    rootsValid = true;

    if (util::checkZero(coeffs[N])) {

        return roots = Polynomial<N - 1>(coeffs).solve();
    }

    switch (N) {

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

template <size_t N>
auto m::Polynomial<N>::getCoeff(size_t index) const {

    if (index > N) return comp::fromCartesian(0, 0);

    return coeffs[index];
}

template <size_t N>
void m::Polynomial<N>::setCoeff(size_t index, const m::comp &value) {

    if (index > N) throw std::invalid_argument("m::exception: cannot extend polynomial");

    coeffs[index] = value;
    rootsValid = false;
}

template <size_t N>
auto m::operator+(const m::Polynomial<N> &lhs, const m::comp &rhs) {

    auto result = lhs;

    result.coeffs[0] += rhs;

    return result;
}

template <size_t N>
auto m::operator+(const m::comp &lhs, const m::Polynomial<N> &rhs) {

    return rhs + lhs;
}

template <size_t N, size_t M>
auto m::operator+(const m::Polynomial<N> &lhs, const m::Polynomial<M> &rhs) {

    if (N > M) {

        auto result = lhs;

        for (size_t i = 0; i < M; i++) {

            result.coeffs[i] += rhs.coeffs[i];
        }

        return result;

    } else {

        auto result = rhs;

        for (size_t i = 0; i < N; i++) {

            result.coeffs[i] += lhs.coeffs[i];
        }

        return result;
    }
}

template <size_t N>
auto m::operator-(const m::Polynomial<N> &rhs) {

    auto result = rhs;

    for (size_t i = 0; i < N; i++) {

        result.coeffs[i] = -result.coeffs[i];
    }

    return result;
}

template <size_t N>
auto m::operator-(const m::Polynomial<N> &lhs, const m::comp &rhs) {

    return lhs + (-rhs);
}

template <size_t N>
auto m::operator-(const m::comp &lhs, const m::Polynomial<N> &rhs) {

    return lhs + (-rhs);
}

template <size_t N, size_t M>
auto m::operator-(const m::Polynomial<N> &lhs, const m::Polynomial<M> &rhs) {

    return lhs + (-rhs);
}

template <size_t N>
auto m::operator*(const m::Polynomial<N> &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < N; i++) {

        result.coeffs[i] *= rhs;
    }

    return result;
}

template <size_t N>
auto m::operator*(const m::comp &lhs, const m::Polynomial<N> &rhs) {

    return rhs * lhs;
}

template <size_t N, size_t M>
auto m::operator*(const m::Polynomial<N> &lhs, const m::Polynomial<M> &rhs) {

    Polynomial<N + M> result;

    for (size_t i = 0; i < N + M; i++) {

        comp c = 0;

        for (size_t j = 0; j <= i; j++) {

            size_t k = i - j;

            auto a = j > N ? 0 : lhs.coeffs[j];
            auto b = k > M ? 0 : rhs.coeffs[k];

            c += a * b;
        }

        result.coeffs[i] = c;
    }

    return result;
}

template <size_t N>
auto m::operator/(const m::Polynomial<N> &lhs, const m::comp &rhs) {

    auto result = lhs;

    for (size_t i = 0; i < N; i++) {

        result.coeffs[i] /= rhs;
    }

    return result;
}

template <size_t N, size_t M>
auto m::operator==(const Polynomial<N> &lhs, const Polynomial<M> &rhs) {

    auto lDeg = lhs.actualDegree();
    auto rDeg = rhs.actualDegree();

    if (!util::checkEqual(lDeg, rDeg)) return false;

    for (size_t i = 0; i < lDeg; i++) {

        if (!util::checkEqual(lhs.coeffs[i], rhs.coeffs[i])) return false;
    }

    return true;
}

template <size_t N, size_t M>
auto m::operator!=(const Polynomial<N> &lhs, const Polynomial<M> &rhs) {


    return !(lhs == rhs);
}

template <size_t N>
auto &m::operator<<(std::ostream &lhs, const m::Polynomial<N> &rhs) {

    lhs << std::fixed << std::setprecision(m_PRECISION);

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

    if (!nonZero) lhs << "0";

    return lhs;
}

template <size_t N>
m::Polynomial<N - 1> m::differentiate(const m::Polynomial<N> &polynomial) {

    Polynomial<N - 1> result;

    for (size_t i = 1; i < N + 1; i++) {

        result.setCoeff(i - 1, m::comp::fromCartesian(i, 0) * polynomial.getCoeff(i));
    }

    return result;
}

template <size_t N>
m::Polynomial<N + 1> m::integrate(const m::Polynomial<N> &polynomial) {

    Polynomial<N + 1> result;
    result.setCoeff(0, m::comp::fromCartesian(0, 0));

    for (size_t i = 0; i < N + 1; i++) {

        result.setCoeff(i + 1, polynomial.getCoeff(i) / m::comp::fromCartesian(i + 1, 0)); 
    }

    return result;
}

#define __m_impl__
#endif
