
template <typename ...Q>
mth::ComplexSolutions mth::ComplexSolutions::finite(Q... args) noexcept {

    return ComplexSolutions(std::unordered_set<mth::comp>({args...}));
}

template <typename ...Q>
mth::Polynomial mth::Polynomial::fromCoeffs(Q... args) {
    return fromCoeffs(std::vector<comp> {static_cast<comp>(args)...});
}
