
#ifdef __mth_mat_content_toggle__

#if defined(N) && defined(M)

template <typename T>
class tmat<T, N, M> {

#elif defined(N)

template <typename T, size_t M>
class tmat<T, N, M> {

#elif defined(M)

template <typename T, size_t N>
class tmat<T, N, M> {

#else

template <typename T, size_t N, size_t M>
class tmat {

#endif

private:

    std::array<T, N * M> values;

    // 0,0 is top left
    constexpr static size_t getIndex(size_t x, size_t y) noexcept {

#ifdef mth_ROW_MAJOR

        return x + y * N;

#else

        return x * M + y;   

#endif

    }

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
    constexpr tmat(const std::array<tvec<T, N>, M> &rows) noexcept {

        for (size_t i = 0; i < M; i++) {

            setRow(i, rows[i]);
        }
    }

    constexpr tmat(const std::array<T, N * M> &values) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                // Always interpret row-major
                get(x, y) = values[x + y * N];
            }
        }
    }

    // Initialize with M rows as tvec<T, N>
    template <typename ...Q, typename std::enable_if<sizeof...(Q) == M, int>::type = 0>
    constexpr tmat(Q... args) noexcept 
        :tmat(std::array<tvec<T, N>, M>{static_cast<tvec<T, N>>(args)...}) {}

#endif

    // Initialize with N*M values written out row-major so that it looks nice with line breaks and hides implementation
    template <typename ...Q, typename std::enable_if<sizeof...(Q) == N * M, int>::type = 0>
    constexpr tmat(Q... args) noexcept
        :values{static_cast<T>(args)...} {

#ifndef mth_ROW_MAJOR

        this->values = static_cast<std::array<T, N * M>>(transpose().values);

#endif

    }

    template <typename U>
    constexpr operator tmat<U, N, M>() const noexcept {

        tmat<U, N, M> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                result.get(x, y) = static_cast<U>(get(x, y));
            }
        }

        return result;
    }

    constexpr size_t size() const noexcept {

        return N * M;
    }

// Convenient getters for tmat<T, 1, 1> without unnecessary indicces
// get() then converts a 1x1 matrix to a scalar
#if defined(N) && defined(M)

    constexpr const T &get() const noexcept {

        return get(0, 0);
    }

    constexpr T &get() noexcept {

        return get(0, 0);
    }

#endif

    constexpr const T &get(size_t x, size_t y) const noexcept {

        return values[getIndex(x, y)];
    }

    constexpr T &get(size_t x, size_t y) noexcept {

        return values[getIndex(x, y)];
    }

    // TODO: Look into getting references here
    constexpr tvec<T, N> getRow(size_t y) const noexcept {

        tvec<T, N> result;

        for (size_t i = 0; i < N; i++) {

            result.get(i) = get(i, y);
        }

        return result;
    }

    constexpr tvec<T, M> getColumn(size_t x) const noexcept {

        tvec<T, M> result;

        for (size_t i = 0; i < M; i++) {
            
            result.get(i) = get(x, i);
        }

        return result;
    }

    constexpr void setRow(size_t y, const tvec<T, N> &value) noexcept {

        for (size_t x = 0; x < N; x++) {

            get(x, y) = value.get(x);
        }
    }

    constexpr void setColumn(size_t x, const tvec<T, M> &value) noexcept {

        for (size_t y = 0; y < M; y++) {

            get(x, y) = value.get(y);
        }
    }

    constexpr std::array<tvec<T, N>, M> rows() const noexcept {

        std::array<tvec<T, N>, M> result;

        for (size_t y = 0; y < M; y++) {

            result[y] = getRow(y);
        }

        return result;
    }

    constexpr std::array<tvec<T, M>, N> columns() const noexcept {

        std::array<tvec<T, M>, N> result;

        for (size_t x = 0; x < N; x++) {

            result[x] = getColumn(x);
        }

        return result;
    }

// Minor matrices only exist for matrices with multiple rows and multiple columns
#if !defined(N) && !defined(M)

    // Returns the matrix found by removing the specified column and row
    constexpr tmat<T, N - 1, M - 1> minor(size_t x, size_t y) const noexcept {

        tmat<T, N - 1, M - 1> result;

        size_t ry = 0;

        for (size_t iy = 0; iy < M; iy++) {

            if (iy == y) continue;

            size_t rx = 0;

            for (size_t ix = 0; ix < N; ix++) {

                if (ix == x) continue;

                result.get(rx, ry) = get(ix, iy);

                rx++;
            }

            ry++;
        }

        return result;
    }

#endif

    // TODO: Update tmat_aug to store determinant multipliers so that echelon form can be used for this
    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr T det() const noexcept {

#if defined(N) || defined(M)

        return get(0, 0);

#else

        T result = 0;
        T s = 1;

        // Can be any row
        size_t row = 0;

        for (size_t x = 0; x < N; x++) {

            auto c = minor(x, row).det();
            result += s * get(x, row) * c;
            s *= -1;
        }

        return result;

#endif

    }

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr bool singular() const noexcept {

        return util::isZero(det());
    }

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr tmat<T, N, M> cofactors() const noexcept {

#if defined(N) || defined(M)

        return tmat<T, 1, 1>::identity();

#else

        tmat<T, N, N> result;
        T s = 1;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                auto c = minor(x, y).det();
                result.get(x, y) = s * c;
                s *= -1;
            }

            if (N % 2 == 0) s *= -1;
        }

        return result;

#endif

    }

    constexpr tmat<T, M, N> transpose() const noexcept {

        tmat<T, M, N> result;

        for (size_t x = 0; x < M; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = get(y, x);
            }
        }

        return result;
    }

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr tmat<T, N, M> adjoint() const noexcept {

        return cofactors().transpose();
    }

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr tmat<T, N, M> inverse() const noexcept {

#ifdef mth_ELIMINATION

        auto id = identity().rows();

        tmat_aug<T, N, tvec<T, N>> augmented(*this, id);

        return tmat<T, N, N>(augmented.solve());

#else

        auto determinant = det();

        return adjoint() / determinant;

#endif 

    }

    // Returns *this / det()
    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr tmat<T, N, M> unit() const noexcept {

        auto determinant = det();

        return *this / determinant;
    }

    template <size_t n = N, size_t m = M, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type = 0>
    constexpr static tmat<T, N, M> identity() noexcept {
        
        tmat<T, N, N> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = x == y ? 1 : 0;
            }
        }

        return result;
    }

    constexpr tmat<T, N, M> &operator+=(const tmat<T, N, M> &rhs) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) += rhs.get(x, y);
            }
        }

        return *this;
    }

    constexpr tmat<T, N, M> &operator-=(const tmat<T, N, M> &rhs) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) -= rhs.get(x, y);
            }
        }

        return *this;
    }

    constexpr tmat<T, N, M> &operator*=(const T &rhs) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) *= rhs;
            }
        }

        return *this;
    }

    constexpr tmat<T, N, M> &operator/=(const T &rhs) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                this->get(x, y) /= rhs;
            }
        }

        return *this;
    }

};

#endif

