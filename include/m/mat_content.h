
#ifdef __m_mat_content_toggle__

// template <typename T, size_t N, size_t M>
// class tmat {

private:

    std::array<T, N * M> values;

    // 0,0 is top left
    static size_t getIndex(size_t x, size_t y);

public:

    // Friend other template instantiations
    template <typename _T, size_t _N, size_t _M>
    friend class tmat;

    // Default initializes to zero
    constexpr tmat() noexcept {

        values.fill(0);
    }

    // These are ambiguous / redundant for 1x1 matrices
#if !defined(N) || !defined(M)

    // Constructors take a list of rows, never a list of columns for consistency with value lists
    tmat(const std::array<tvec<T, N>, M> &rows) noexcept;

    tmat(const std::array<T, N * M> &values) noexcept;

    // Initialize with M rows as tvec<T, N>
    template <typename ...Q, typename std::enable_if<sizeof...(Q) == M, int>::type = 0>
    tmat(Q... args) noexcept 
        :tmat(std::array<tvec<T, N>, M>{static_cast<tvec<T, N>>(args)...}) {}

#endif

    // Initialize with N*M values written out row-major so that it looks nice with line breaks and hides implementation
    template <typename ...Q, typename std::enable_if<sizeof...(Q) == N * M, int>::type = 0>
    constexpr tmat(Q... args) noexcept
        :values{static_cast<T>(args)...} {

#ifndef m_ROW_MAJOR

        this->values = static_cast<std::array<T, N * M>>(transpose().values);

#endif

    }

    constexpr size_t size() const noexcept {

        return N * M;
    }

// Convenient getters for tmat<T, 1, 1> without unnecessary indicces
// get() then converts a 1x1 matrix to a scalar
#if defined(N) && defined(M)

    const T &get() const;
    T &get();

#endif

    const T &get(size_t x, size_t y) const;

    T &get(size_t x, size_t y);

    // TODO: Look into getting references here
    tvec<T, N> getRow(size_t y) const;
    tvec<T, M> getColumn(size_t x) const;

    void setRow(size_t y, const tvec<T, N> &value);
    void setColumn(size_t x, const tvec<T, M> &value);

    std::array<tvec<T, N>, M> rows() const;
    std::array<tvec<T, M>, N> columns() const;

// Minor matrices only exist for matrices with multiple rows and multiple columns
#if !defined(N) && !defined(M)

    // Returns the matrix found by removing the specified column and row
    tmat<T, N - 1, M - 1> minor(size_t x, size_t y) const;

#endif

    // TODO: Update tmat_aug to store determinant multipliers so that echelon form can be used for this
    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    T det() const noexcept;

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    bool singular() const noexcept;

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    tmat<T, N, M> cofactors() const;

    tmat<T, M, N> transpose() const;

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    tmat<T, N, M> adjoint() const;

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    tmat<T, N, M> inverse() const;

    // Returns *this / det()
    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    tmat<T, N, M> unit() const;

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    static tmat<T, N, M> identity();

    tmat<T, N, M> &operator+=(const tmat<T, N, M> &rhs);
    tmat<T, N, M> &operator-=(const tmat<T, N, M> &rhs);
    tmat<T, N, M> &operator*=(const T &rhs);
    tmat<T, N, M> &operator/=(const T &rhs);
// };

#endif

