
template <typename ...Q>
m::ComplexSolutions m::ComplexSolutions::finite(Q... args) noexcept {

    return ComplexSolutions(std::unordered_set<m::comp>({args...}));
}

template <typename ...Q>
m::Polynomial m::Polynomial::fromCoeffs(Q... args) {
    return fromCoeffs(std::vector<comp> {static_cast<comp>(args)...});
}
