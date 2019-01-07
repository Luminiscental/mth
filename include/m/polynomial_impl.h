
template <typename ...Q>
m::ComplexSolutions m::ComplexSolutions::finite(Q... args) noexcept {

    return ComplexSolutions(std::unordered_set<m::comp>({args...}));
}

template <typename ...Q>
m::Polynomial::Polynomial(Q... args) noexcept
    :Polynomial(std::vector<comp> {static_cast<comp>(args)...}) {}
