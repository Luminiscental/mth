
#ifdef __Maths_mat_content_toggle__

// template <typename T, size_t N>
// class tmat {

private:
    std::array<T, N * N> values;

    static size_t getIndex(size_t x, size_t y) { // NOTE: 0,0 is top left

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > N - 1) throw std::out_of_range("m::exception: row index out of bounds");

#ifdef Maths_ROW_MAJOR // TODO: Maybe make this choosable per-object or have a different class for row-major matrices
        return x + y * N;
#else
        return x * N + y;   
#endif
    }

public:
    const T &get(size_t x, size_t y) const {

        return values[getIndex(x, y)];
    }

    T &get(size_t x, size_t y) {

        return values[getIndex(x, y)];
    }

    tvec<T, N> getRow(size_t y) const {

        tvec<T, N> result;

        for (size_t i = 0; i < N; i++) {

            result.get(i) = get(i, y);
        }

        return result;
    }

    tvec<T, N> getColumn(size_t x) const {

        tvec<T, N> result;

        for (size_t i = 0; i < N; i++) {
            
            result.get(i) = get(x, i);
        }

        return result;
    }

    tmat() {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                get(x, y) = (x == y) ? 1 : 0;
            }
        }
    }

#ifdef Maths_ROW_MAJOR

    tmat(std::array<tvec<T, N>, N> rows) {

        for (size_t y = 0; y < N; y++) {

            for (size_t x = 0; x < N: x++) {

                get(x, y) = rows[y].get(x);
            }
        }
    }

#else

    tmat(std::array<tvec<T, N>, N> columns) {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                get(x, y) = columns[x].get(y);
            }
        }
    }

#endif

    tmat(const std::array<T, N * N> &values) : values(values) {}

    tmat(const tmat<T, N> &other) : tmat(other.values) {}

    template <typename ...Q, typename std::enable_if<sizeof...(Q) == N * N, int>::type = 0>
    tmat(Q... args) : values{static_cast<T>(args)...} {}

#ifndef __Maths_mat_basecaseimpl__
    tmat<T, N - 1> minor(size_t x, size_t y) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > N - 1) throw std::out_of_range("m::exception: row index out of bounds");

        tmat<T, N - 1> result;

        size_t ry = 0;

        for (size_t iy = 0; iy < N; iy++) {

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
#else
    T minor(size_t x, size_t y) const {

        if (x > 1) throw std::out_of_range("m::exception: column index out of bounds");
        if (y > 1) throw std::out_of_range("m::exception: row index out of bounds");

        return get((x + 1) % 2, (y + 1) % 2);
    }
#endif

    T det() const {

        T result = 0;
        T s = 1;
        size_t row = 0; // can be any row

        for (size_t x = 0; x < N; x++) {

#ifndef __Maths_mat_basecaseimpl__
            T c = minor(x, row).det();
#else
            T c = minor(x, row);
#endif
            result += s * get(x, row) * c;
            s *= -1;
        }

        return result;
    }

    tmat<T, N> cofactors() const {

        tmat<T, N> result;
        T s = 1;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

#ifndef __Maths_mat_basecaseimpl__
                T c = minor(x, y).det();
#else
                T c = minor(x, y);
#endif
                result.get(x, y) = s * c;
                s *= -1;
            }

            if (N % 2 == 0) s *= -1;
        }

        return result;
    }

    tmat<T, N> transpose() const {

        tmat<T, N> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = get(y, x);
            }
        }

        return result;
    }

    tmat<T, N> adjoint() const {

        return cofactors().transpose();
    }

    tmat<T, N> inverse() const {

        T determinant = det();

        if (util::checkZero(determinant)) throw std::invalid_argument("m::exception: inverse() called on singular matrix");

        return adjoint() / determinant;
    }

    tmat<T, N> unit() const {

        T determinant = det();

        if (util::checkZero(determinant)) throw std::invalid_argument("m::exception: unit() called on singular matrix");

        return *this / determinant;
    }

    static tmat<T, N> identity() {
        
        return tmat<T, N>(); 
    }

    friend tmat<T, N> operator+(const tmat<T, N> &a, const tmat<T, N> &b) {

        tmat<T, N> result(a);

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) += b.get(x, y);
            }
        }

        return result;
    }

    friend tmat<T, N> operator-(const tmat<T, N> &a, const tmat<T, N> &b) {

        tmat<T, N> result(a);

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) -= b.get(x, y);
            }
        }

        return result;
    }

    friend tmat<T, N> operator*(T scalar, const tmat<T, N> &matrix) {

        tmat<T, N> result(matrix);

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) *= scalar;
            }
        }

        return result;
    }

    friend tmat<T, N> operator*(const tmat<T, N> &matrix, T scalar) {

        return scalar * matrix;
    }

    friend tmat<T, N> operator/(const tmat<T, N> &matrix, T scalar) {

        tmat<T, N> result(matrix);

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) /= scalar;
            }
        }

        return result;
    }

    friend tmat<T, N> operator*(const tmat<T, N> &a, const tmat<T, N> &b) {

        tmat<T, N> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(x, y) = a.getRow(y).dot(b.getColumn(x));
            }
        }

        return result;
    }

    friend tmat<T, N> operator/(const tmat<T, N> &a, const tmat<T, N> &b) {

        return b.inverse() * a;
    }

    friend tvec<T, N> operator*(const tmat<T, N> &matrix, const tvec<T, N> &vector) {

        tvec<T, N> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < N; y++) {

                result.get(y) = matrix.getRow(y).dot(vector);
            }
        }

        return result;
    }

    // TODO: Either make this less crappy or remove it
    friend std::ostream &operator<<(std::ostream &stream, const tmat<T, N> &matrix) {

        stream << std::fixed << std::setprecision(2);

        stream << std::endl << "__";
        
        for (size_t i = 0; i < N + 1; i++) {

            stream << "\t";
        }

        stream << "\b \b__" << std::endl;

        for (size_t y = 0; y < N; y++) {

            stream << "|\t";

            for (size_t x = 0; x < N - 1; x++) {

                stream << matrix.get(x, y) << "\t";
            }

            stream << matrix.get(N - 1, y) << "\t|" << std::endl;
        }

        stream << "¯¯";

        for (size_t i = 0; i < N + 1; i++) {

            stream << "\t";
        }

        return stream << "\b \b¯¯";
    }
// };

#endif

