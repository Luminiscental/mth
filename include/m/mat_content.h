
#ifdef __m_mat_content_toggle__

// template <typename T, size_t N, size_t M>
// class tmat {

private:

    std::array<T, N * M> values;

    static auto getIndex(size_t x, size_t y) { // NOTE: 0,0 is top left

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

#ifdef m_ROW_MAJOR

        return x + y * N;

#else

        return x * M + y;   

#endif

    }

public:

    constexpr tmat() noexcept {

        values.fill(0);
    }

    // NOTE: Constructors take a list of rows, never a list of columns for consistency with value lists
    tmat(const std::array<tvec<T, N>, M> &rows) noexcept {

        for (size_t i = 0; i < M; i++) {

            setRow(i, rows[i]);
        }
    }

    tmat(const std::array<T, N * M> &values) noexcept
        :values(values) {

#ifndef m_ROW_MAJOR

        this->values = transpose().values; // NOTE: Value lists are row major

#endif

    }

    template <typename ...Q, typename std::enable_if<sizeof...(Q) == M, int>::type = 0>
    tmat(Q... args) noexcept {

        std::array<tvec<T, N>, M> rows{static_cast<tvec<T, N>>(args)...};

        for (size_t i = 0; i < M; i++) {

            setRow(i, rows[i]);
        }
    }

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

    const auto &get(size_t x, size_t y) const {

        return values[getIndex(x, y)];
    }

    auto &get(size_t x, size_t y) {

        return values[getIndex(x, y)];
    }

    // TODO: Look into getting references here
    auto getRow(size_t y) const {

        tvec<T, N> result;

        for (size_t i = 0; i < N; i++) {

            result.get(i) = get(i, y);
        }

        return result;
    }

    auto getColumn(size_t x) const {

        tvec<T, M> result;

        for (size_t i = 0; i < M; i++) {
            
            result.get(i) = get(x, i);
        }

        return result;
    }

    void setRow(size_t y, const tvec<T, N> &value) {

        for (size_t x = 0; x < N; x++) {

            get(x, y) = value.get(x);
        }
    }

    void setColumn(size_t x, const tvec<T, M> &value) {

        for (size_t y = 0; y < M; y++) {

            get(x, y) = value.get(y);
        }
    }

    auto rows() const {

        std::array<tvec<T, N>, M> result;

        for (size_t y = 0; y < M; y++) {

            result[y] = getRow(y);
        }

        return result;
    }

    auto columns() const {

        std::array<tvec<T, M>, N> result;

        for (size_t x = 0; x < N; x++) {

            result[x] = getColumn(x);
        }

        return result;
    }

#if defined(N)

#if defined(M) // tmat<2, 2>

    auto minor(size_t x, size_t y) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

        return get((x + 1) % 2, (y + 1) % 2);
    }

#else // tmat<2, M>

    auto minor(size_t x, size_t y) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

        return getColumn((x + 1) % 2);
    }

#endif

#elif defined(M) // tmat<N, 2>

    auto minor(size_t x, size_t y) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

        return getRow((y + 1) % 2);
    }

#else // tmat<N, M>

    auto minor(size_t x, size_t y) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

        tmat<T, N - 1, M - 1> result;

        size_t ry = 0;

        for (size_t iy = 0; iy < M; iy++) {

            if (iy != y) {

                size_t rx = 0;

                for (size_t ix = 0; ix < N; ix++) {

                    if (ix != x) {

                        result.get(rx, ry) = get(ix, iy);

                        rx++;
                    }
                }

                ry++;
            }
        }

        return result;
    }

#endif

    // TODO: Update tmat_aug to store determinant multipliers so that echelon form can be used for this
    template <typename std::enable_if<N == M, int>::type = 0>
    auto det() const noexcept {

        T result = 0;
        T s = 1;
        size_t row = 0; // NOTE: can be any row

        for (size_t x = 0; x < N; x++) {

            auto c = minor(x, row)

#ifndef N // enable_if means either tmat<2, 2> or tmat<N, N> 

            .det()

#endif

            ;

            result += s * get(x, row) * c;
            s *= -1;
        }

        return result;
    }

    template <typename std::enable_if<N ==M, int>::type = 0>
    bool singular() const noexcept {

        return util::checkZero(det());
    }

    template <typename std::enable_if<N == M, int>::type = 0>
    auto cofactors() const {

        tmat<T, N, N> result;
        T s = 1;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                auto c = minor(x, y)

#ifndef N // enable_if means either tmat<2, 2> or tmat<N, N>

                .det()

#endif

                ;

                result.get(x, y) = s * c;
                s *= -1;
            }

            if (N % 2 == 0) s *= -1;
        }

        return result;
    }

    auto transpose() const {

        tmat<T, M, N> result;

        for (size_t x = 0; x < M; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = get(y, x);
            }
        }

        return result;
    }

    template <typename std::enable_if<N == M, int>::type = 0>
    auto adjoint() const {

        return cofactors().transpose();
    }

    template <typename std::enable_if<N == M, int>::type = 0>
    auto inverse() const {

#ifdef m_ELIMINATION

        if (singular()) throw std::invalid_argument("m::exception: inverse() called on singular matrix");

        auto id = identity().rows();

        tmat_aug<T, N, tvec<T, N>> augmented(*this, id);

        return tmat<T, N, N>(augmented.reducedRowEchelon().auxilary());

#else

        auto determinant = det();

        if (util::checkZero(determinant)) throw std::invalid_argument("m::exception: inverse() called on singular matrix");

        return adjoint() / determinant;

#endif 

    }

    template <typename std::enable_if<N == M, int>::type = 0>
    auto unit() const {

        auto determinant = det();

        if (util::checkZero(determinant)) throw std::invalid_argument("m::exception: unit() called on singular matrix");

        return *this / determinant;
    }

    template <typename std::enable_if<N == M, int>::type = 0>
    static auto identity() {
        
        tmat<T, N, N> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = x == y ? 1 : 0;
            }
        }

        return result;
    }

    auto &operator+=(const tmat<T, N, M> &rhs) {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) += rhs.get(x, y);
            }
        }

        return *this;
    }

    auto &operator-=(const tmat<T, N, M> &rhs) {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) -= rhs.get(x, y);
            }
        }

        return *this;
    }

    auto &operator*=(const T &rhs) {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) *= rhs;
            }
        }

        return *this;
    }

    auto &operator/=(const T &rhs) {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) /= rhs;
            }
        }

        return *this;
    }

    friend auto operator+(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) {

        auto result = lhs;

        return result += rhs;
    }

    friend auto operator-(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) {

        auto result = lhs;

        return result -= rhs;
    }

    friend auto operator*(const T &lhs, const tmat<T, N, M> &rhs) {

        auto result = rhs;

        return result *= lhs;
    }

    friend auto operator*(const tmat<T, N, M> &lhs, const T &rhs) {

        return rhs * lhs;
    }

    friend auto operator/(const tmat<T, N, M> &lhs, const T &rhs) {

        auto result = lhs;

        return result /= rhs;
    }

    template <size_t O>
    friend auto operator*(const tmat<T, N, M> &lhs, const tmat<T, O, N> &rhs) {

        tmat<T, O, M> result;

        for (size_t x = 0; x < O; x++) {

            for (size_t y = 0; y < M; y++) {

                result.get(x, y) = lhs.getRow(y).dot(rhs.getColumn(x));
            }
        }

        return result;
    }

    template <typename std::enable_if<N == M, int>::type = 0>
    friend auto operator/(const tmat<T, N, N> &lhs, const tmat<T, N, N> &rhs) {

        return lhs * rhs.inverse();
    }

    friend auto operator*(const tmat<T, N, M> &lhs, const tvec<T, N> &rhs) {

        tvec<T, M> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                result.get(y) = vec::dot(lhs.getRow(y), rhs);
            }
        }

        return result;
    }

    // TODO: Either make this less crappy or remove it
    friend auto &operator<<(std::ostream &lhs, const tmat<T, N, M> &rhs) {

        lhs << std::fixed << std::setprecision(m_PRECISION);

        for (size_t y = 0; y < M; y++) {

            lhs << "|\t";

            for (size_t x = 0; x < N - 1; x++) {

                lhs << rhs.get(x, y) << "\t";
            }

            lhs << rhs.get(N - 1, y) << "\t|";

            if (y < M - 1) lhs << std::endl;
        }

        return lhs;
    }

// };

#endif

