
#ifdef __m_mat_content_toggle__

// template <typename T, size_t N, size_t M>
// class tmat {

private:

    std::array<T, N * M> values;

    static auto getIndex(size_t x, size_t y); // NOTE: 0,0 is top left

public:

    constexpr tmat() noexcept {

        values.fill(0);
    }

    // NOTE: Constructors take a list of rows, never a list of columns for consistency with value lists
    tmat(const std::array<tvec<T, N>, M> &rows) noexcept;

    tmat(const std::array<T, N * M> &values) noexcept;

    template <typename ...Q, typename std::enable_if<sizeof...(Q) == M, int>::type = 0>
    tmat(Q... args) noexcept 
        :tmat(std::array<tvec<T, N>, M>{static_cast<tvec<T, N>>(args)...}) {}

    template <typename ...Q, typename std::enable_if<sizeof...(Q) == N * M, int>::type = 0>
    constexpr tmat(Q... args) noexcept
        :values{static_cast<T>(args)...} {

#ifndef m_ROW_MAJOR

        this->values = transpose().values;

#endif

    }

    constexpr size_t size() const noexcept {

        return N * M;
    }

    const T &get(size_t x, size_t y) const;

    T &get(size_t x, size_t y);

    // TODO: Look into getting references here
    auto getRow(size_t y) const;

    auto getColumn(size_t x) const;

    void setRow(size_t y, const tvec<T, N> &value);

    void setColumn(size_t x, const tvec<T, M> &value);

    auto rows() const;

    auto columns() const;

#if defined(N)

#if defined(M) // 2,2

    T minor(size_t x, size_t y) const;

#else // 2,M

    tvec<T, M - 1> minor(size_t x, size_t y) const;

#endif

#else

#if defined(M) // N,2

    tvec<T, N - 1> minor(size_t x, size_t y) const;

#else // N,M

    tmat<T, N - 1, M - 1> minor(size_t x, size_t y) const;

#endif

#endif

    // TODO: Update tmat_aug to store determinant multipliers so that echelon form can be used for this
    template <typename std::enable_if<N == M, int>::type = 0>
    auto det() const noexcept;

    template <typename std::enable_if<N ==M, int>::type = 0>
    bool singular() const noexcept;

    template <typename std::enable_if<N == M, int>::type = 0>
    auto cofactors() const;

    auto transpose() const;

    template <typename std::enable_if<N == M, int>::type = 0>
    auto adjoint() const;

    template <typename std::enable_if<N == M, int>::type = 0>
    auto inverse() const;

    template <typename std::enable_if<N == M, int>::type = 0>
    auto unit() const;

    template <typename std::enable_if<N == M, int>::type = 0>
    static auto identity();

    auto &operator+=(const tmat<T, N, M> &rhs);

    auto &operator-=(const tmat<T, N, M> &rhs);

    auto &operator*=(const T &rhs);

    auto &operator/=(const T &rhs);
// };

#endif

