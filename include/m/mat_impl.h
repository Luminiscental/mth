
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>

// tmat_aug is not specialized so only include it in the general case
#if !defined(N) && !defined(M)

template <typename T, size_t N, typename A>
m::tmat_aug<T, N, A>::tmat_aug() noexcept {}

template <typename T, size_t N, typename A>
m::tmat_aug<T, N, A>::tmat_aug(const m::tmat<T, N, N> &matrix, const std::array<A, N> &aux) noexcept
    :matrix(matrix), aux(aux) {}

template <typename T, size_t N, typename A>
m::tmat<T, N, N> m::tmat_aug<T, N, A>::coefficients() const {

    return matrix;
}

template <typename T, size_t N, typename A>
m::tvec<A, N> m::tmat_aug<T, N, A>::auxilary() const {

    return aux;
}

template <typename T, size_t N, typename A>
m::tvec<A, N> m::tmat_aug<T, N, A>::solve() const {

    if (singular()) throw std::invalid_argument("m::exception: solve() called on singular system");

    return reducedRowEchelon().auxilary();
}

template <typename T, size_t N, typename A>
size_t m::tmat_aug<T, N, A>::leadingIndex(size_t row) const {

    for (size_t i = 0; i < N; i++) {

        if (!util::isZero(matrix.get(i, row))) return i;
    }

    return N;
}

template <typename T, size_t N, typename A>
T m::tmat_aug<T, N, A>::leadingValue(size_t row) const {

    auto index = leadingIndex(row);

    if (index == N) return 0;

    return matrix.get(index, row);
}

template <typename T, size_t N, typename A>
bool m::tmat_aug<T, N, A>::columnIsZero(size_t x) const {

    if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
    
    for (size_t y = 0; y < N; y++) {

        if (!util::isZero(matrix.get(x, y))) return false;
    }

    return true;
}

template <typename T, size_t N, typename A>
bool m::tmat_aug<T, N, A>::rowIsZero(size_t y) const {

    if (y > N - 1) throw std::out_of_range("m::exception: row index out of bounds");

    for (size_t x = 0; x < N; x++) {

        if (!util::isZero(matrix.get(x, y))) return false;
    }

    return true;
}

template <typename T, size_t N, typename A>
bool m::tmat_aug<T, N, A>::hasZeroRow() const {

    for (size_t y = 0; y < N; y++) {

        if (rowIsZero(y)) return true;
    }

    return false;
}

template <typename T, size_t N, typename A>
bool m::tmat_aug<T, N, A>::singular() const {

    return rowEchelon().hasZeroRow();
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::swapRows(size_t a, size_t b) {

    auto rowA = matrix.getRow(a);
    auto rowB = matrix.getRow(b);

    matrix.setRow(a, rowB);
    matrix.setRow(b, rowA);

    auto auxA = aux[a];
    auto auxB = aux[b];

    aux[a] = auxB;
    aux[b] = auxA;
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::scaleRow(size_t index, T scalar) {

    matrix.setRow(index, scalar * matrix.getRow(index));
    aux[index] = scalar * aux[index];
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::addRow(size_t targetRow, size_t sourceRow, T scalar) {

    auto addValue = scalar * matrix.getRow(sourceRow);
    auto resultValue = matrix.getRow(targetRow) + addValue;

    matrix.setRow(targetRow, resultValue);

    aux[targetRow] += scalar * aux[sourceRow];
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::setRow(size_t index, const m::tvec<T, N> &val, A auxVal) {

    matrix.setRow(index, val);
    aux[index] = auxVal;
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::eliminateFromBelow(size_t x, size_t y) {

    // Want to add a row scaled to have this value in column x
    auto targetValue = -matrix.get(x, y);

    if (util::isZero(targetValue)) return;

    // Look for a row to use below y
    for (size_t iy = y + 1; iy < N; iy++) {

        // Continue if values to the left would be affected
        if (leadingIndex(iy) < x) continue;

        auto val = matrix.get(x, iy);

        // If it is zero it can't become targetValue
        if (!util::isZero(val)) {

            addRow(y, iy, targetValue / val); 
            return;
        }
    }

    throw std::invalid_argument("m::exception: eliminateFromBelow() called on non-eliminable element");
}

template <typename T, size_t N, typename A>
void m::tmat_aug<T, N, A>::eliminateFromRight(size_t x, size_t y) {

    // Want to add a row scaled to have this value in column x
    auto targetValue = -matrix.get(x, y);

    if (util::isZero(targetValue)) return;

    for (size_t iy = 0; iy < N; iy++) {

        // We don't want to remove the row
        if (iy == y) continue;

        auto value = matrix.get(x, iy);

        // If the value is zero it can't be made into targetValue
        // Can't affect values to the left of x
        if (!util::isZero(value) && leadingIndex(iy) >= x) {

            addRow(y, iy, targetValue / value);
            return;
        }
    }

    throw std::invalid_argument("m::exception: eliminateFromRight() called on non-eliminable element");
}

template <typename T, size_t N, typename A>
m::tmat_aug<T, N, A> m::tmat_aug<T, N, A>::ordered() const {

    using std::begin;
    using std::end;

    // Store the current ordering of rows, i.e. 0, 1, 2, 3, ...
    std::array<size_t, N> rowIndices;
    std::iota(begin(rowIndices), end(rowIndices), 0);

    // Permute the rowIndices to sort by leadingIndex
    auto compareLeadingIndex = [=] (auto a, auto b) -> bool { return leadingIndex(a) < leadingIndex(b); };
    std::sort(begin(rowIndices), end(rowIndices), compareLeadingIndex);

    tmat_aug<T, N, A> result{matrix, aux};

    // Apply the permutation
    for (size_t i = 0; i < N; i++) {

        result.setRow(i, matrix.getRow(rowIndices[i]), aux[rowIndices[i]]);
    }

    return result;
}

template <typename T, size_t N, typename A>
m::tmat_aug<T, N, A> m::tmat_aug<T, N, A>::rowEchelon() const {

    auto result = ordered();

    // Iterate over columns and eliminate below the diagonal
    for (size_t x = 0; x < N - 1; x++) { 

        // Leave empty columns
        if (result.columnIsZero(x)) continue;

        // Iterate over each row and eliminate the value in that row
        for (size_t y = x + 1; y < N; y++) {

            // If the row is zero we are already done
            if (result.rowIsZero(y)) break;

            // If the value is zero we are already done
            if (util::isZero(result.matrix.get(x, y))) continue;

            result.eliminateFromRight(x, y);

            // If we have eliminated re-order 
            result = result.ordered();

            // Go back to the first row so none are missed
            y = x;
        }
    }

    return result.ordered(); 
}

template <typename T, size_t N, typename A>
m::tmat_aug<T, N, A> m::tmat_aug<T, N, A>::reducedRowEchelon() const {

    // Start from echelon form so that half of the values to eliminate are already zero
    auto result = rowEchelon();

    if (result.hasZeroRow()) throw std::invalid_argument("m::exception: reducedRowEchelon() called on singular matrix");

    for (size_t y = 0; y < N; y++) { 

        auto leadingValue = result.leadingValue(y);

        // Make the diagonals 1
        result.scaleRow(y, 1 / leadingValue);

        // Eliminate values above the diagonal
        for (size_t x = result.leadingIndex(y) + 1; x < N; x++) {

            result.eliminateFromBelow(x, y);
        }
    }

    return result;
}

template <typename T, size_t N, typename A>
std::ostream &m::operator<<(std::ostream &lhs, const m::tmat_aug<T, N, A> &rhs) {

    lhs << std::fixed << std::setprecision(m_PRECISION);

    for (size_t y = 0; y < N; y++) {

        lhs << "|\t";

        for (size_t x = 0; x < N - 1; x++) {

            lhs << rhs.matrix.get(x, y) << "\t";
        }

        lhs << rhs.matrix.get(N - 1, y) << "\t|\t" << rhs.aux[y] << "\t|";

        if (y < N - 1) lhs << std::endl;
    }
    
    return lhs;
}

template <typename T, size_t N, size_t M>
size_t std::hash<m::tmat<T, N, M>>::operator()(const m::tmat<T, N, M> &x) {

    size_t result = 0;

    for (auto row : x.rows()) {

        result ^= hash<decltype(row)>()(row);
    }

    return result;
}

#endif

// Defines for different template specializations

#if defined(N) && defined(M)

#define TEMPLATE_NM template <typename T>
#define TEMPLATE_NN template <typename T> template <size_t n, size_t m, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type>
#define TEMPLATE_NMO template <typename T, size_t O>

#elif defined(N)

#define TEMPLATE_NM template <typename T, size_t M>
#define TEMPLATE_NN template <typename T, size_t M> template <size_t n, size_t m, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type>
#define TEMPLATE_NMO template <typename T, size_t M, size_t O>

#elif defined(M)

#define TEMPLATE_NM template <typename T, size_t N>
#define TEMPLATE_NN template <typename T, size_t N> template <size_t n, size_t m, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type>
#define TEMPLATE_NMO template <typename T, size_t N, size_t O>

#else

#define TEMPLATE_NM template <typename T, size_t N, size_t M>
#define TEMPLATE_NN template <typename T, size_t N, size_t M> template <size_t n, size_t m, typename std::enable_if<(n == m) && (n == N) && (m == M), int>::type>
#define TEMPLATE_NMO template <typename T, size_t N, size_t M, size_t O>

#endif

TEMPLATE_NM
size_t m::tmat<T, N, M>::getIndex(size_t x, size_t y) {

    if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
    if (y > M - 1) throw std::out_of_range("m::exception: row index out of bounds");

#ifdef m_ROW_MAJOR

    return x + y * N;

#else

    return x * M + y;   

#endif

}

#if !defined(N) || !defined(M)

TEMPLATE_NM
m::tmat<T, N, M>::tmat(const std::array<m::tvec<T, N>, M> &rows) noexcept {

    for (size_t i = 0; i < M; i++) {

        setRow(i, rows[i]);
    }
}

TEMPLATE_NM
m::tmat<T, N, M>::tmat(const std::array<T, N * M> &values) noexcept {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            // Always interpret row-major
            get(x, y) = values[x + y * N];
        }
    }
}

#endif

#if defined(N) && defined(M)

TEMPLATE_NM
const T &m::tmat<T, N, M>::get() const {

    return get(0, 0);
}

TEMPLATE_NM
T &m::tmat<T, N, M>::get() {

    return get(0, 0);
}

#endif

TEMPLATE_NM
const T &m::tmat<T, N, M>::get(size_t x, size_t y) const {

    return values[getIndex(x, y)];
}

TEMPLATE_NM
T &m::tmat<T, N, M>::get(size_t x, size_t y) {

    return values[getIndex(x, y)];
}

// TODO: Look into getting references here
TEMPLATE_NM
m::tvec<T, N> m::tmat<T, N, M>::getRow(size_t y) const {

    tvec<T, N> result;

    for (size_t i = 0; i < N; i++) {

        result.get(i) = get(i, y);
    }

    return result;
}

TEMPLATE_NM
m::tvec<T, M> m::tmat<T, N, M>::getColumn(size_t x) const {

    tvec<T, M> result;

    for (size_t i = 0; i < M; i++) {
        
        result.get(i) = get(x, i);
    }

    return result;
}

TEMPLATE_NM
void m::tmat<T, N, M>::setRow(size_t y, const m::tvec<T, N> &value) {

    for (size_t x = 0; x < N; x++) {

        get(x, y) = value.get(x);
    }
}

TEMPLATE_NM
void m::tmat<T, N, M>::setColumn(size_t x, const m::tvec<T, M> &value) {

    for (size_t y = 0; y < M; y++) {

        get(x, y) = value.get(y);
    }
}

TEMPLATE_NM
std::array<m::tvec<T, N>, M> m::tmat<T, N, M>::rows() const {

    std::array<tvec<T, N>, M> result;

    for (size_t y = 0; y < M; y++) {

        result[y] = getRow(y);
    }

    return result;
}

TEMPLATE_NM
std::array<m::tvec<T, M>, N> m::tmat<T, N, M>::columns() const {

    std::array<tvec<T, M>, N> result;

    for (size_t x = 0; x < N; x++) {

        result[x] = getColumn(x);
    }

    return result;
}

#if !defined(N) && !defined(M)

TEMPLATE_NM
m::tmat<T, N - 1, M - 1> m::tmat<T, N, M>::minor(size_t x, size_t y) const {

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

#if defined(N) || defined(M)

TEMPLATE_NN
T m::tmat<T, N, M>::det() const noexcept {

    return get(0, 0);
}

#else

// TODO: Update tmat_aug to store determinant multipliers so that echelon form can be used for this
TEMPLATE_NN
T m::tmat<T, N, M>::det() const noexcept {

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
}

#endif

TEMPLATE_NN
bool m::tmat<T, N, M>::singular() const noexcept {

    return util::isZero(det());
}

#if defined(N) || defined(M)

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::cofactors() const {

    return m::tmat<T, 1, 1>::identity();
}

#else

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::cofactors() const {

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
}

#endif

TEMPLATE_NM
m::tmat<T, M, N> m::tmat<T, N, M>::transpose() const {

    tmat<T, M, N> result;

    for (size_t x = 0; x < M; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = get(y, x);
        }
    }

    return result;
}

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::adjoint() const {

    return cofactors().transpose();
}

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::inverse() const {

#ifdef m_ELIMINATION

    auto id = identity().rows();

    tmat_aug<T, N, tvec<T, N>> augmented(*this, id);

    if (augmented.singular()) throw std::invalid_argument("m::exception: inverse() called on singular matrix");

    return tmat<T, N, N>(augmented.solve());

#else

    auto determinant = det();

    if (util::isZero(determinant)) throw std::invalid_argument("m::exception: inverse() called on singular matrix");

    return adjoint() / determinant;

#endif 

}

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::unit() const {

    auto determinant = det();

    if (util::isZero(determinant)) throw std::invalid_argument("m::exception: unit() called on singular matrix");

    return *this / determinant;
}

TEMPLATE_NN
m::tmat<T, N, M> m::tmat<T, N, M>::identity() {
    
    tmat<T, N, N> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < N; y++) {

            result.get(x, y) = x == y ? 1 : 0;
        }
    }

    return result;
}

TEMPLATE_NM
m::tmat<T, N, M> &m::tmat<T, N, M>::operator+=(const tmat<T, N, M> &rhs) {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            this->get(x, y) += rhs.get(x, y);
        }
    }

    return *this;
}

TEMPLATE_NM
m::tmat<T, N, M> &m::tmat<T, N, M>::operator-=(const tmat<T, N, M> &rhs) {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            this->get(x, y) -= rhs.get(x, y);
        }
    }

    return *this;
}

TEMPLATE_NM
m::tmat<T, N, M> &m::tmat<T, N, M>::operator*=(const T &rhs) {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            this->get(x, y) *= rhs;
        }
    }

    return *this;
}

TEMPLATE_NM
m::tmat<T, N, M> &m::tmat<T, N, M>::operator/=(const T &rhs) {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            this->get(x, y) /= rhs;
        }
    }

    return *this;
}

#undef TEMPLATE_NN
#undef TEMPLATE_NM
#undef TEMPLATE_NMO

#if !defined(N) && !defined(M)

template <typename T, size_t N, size_t M>
m::tmat<T, N, M> m::operator+(const m::tmat<T, N, M> &lhs, const m::tmat<T, N, M> &rhs) {

    auto result = lhs;

    return result += rhs;
}

template <typename T, size_t N, size_t M>
m::tmat<T, N, M> m::operator-(const m::tmat<T, N, M> &lhs, const m::tmat<T, N, M> &rhs) {

    auto result = lhs;

    return result -= rhs;
}

template <typename T, size_t N, size_t M>
m::tmat<T, N, M> m::operator*(const T &lhs, const m::tmat<T, N, M> &rhs) {

    auto result = rhs;

    return result *= lhs;
}

template <typename T, size_t N, size_t M>
m::tmat<T, N, M> m::operator*(const m::tmat<T, N, M> &lhs, const T &rhs) {

    return rhs * lhs;
}

template <typename T, size_t N, size_t M>
m::tmat<T, N, M> m::operator/(const m::tmat<T, N, M> &lhs, const T &rhs) {

    auto result = lhs;

    return result /= rhs;
}

template <typename T, size_t N, size_t M, size_t O>
m::tmat<T, O, M> m::operator*(const m::tmat<T, N, M> &lhs, const m::tmat<T, O, N> &rhs) {

    tmat<T, O, M> result;

    for (size_t x = 0; x < O; x++) {

        for (size_t y = 0; y < M; y++) {

            result.get(x, y) = tvec<T, N>::dot(lhs.getRow(y), rhs.getColumn(x));
        }
    }

    return result;
}

template <typename T, size_t N>
m::tmat<T, N, N> m::operator/(const m::tmat<T, N, N> &lhs, const m::tmat<T, N, N> &rhs) {

    return lhs * rhs.inverse();
}

template <typename T, size_t N, size_t M>
m::tvec<T, M> m::operator*(const m::tmat<T, N, M> &lhs, const m::tvec<T, N> &rhs) {

    tvec<T, M> result;

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            result.get(y) = tvec<T, N>::dot(lhs.getRow(y), rhs);
        }
    }

    return result;
}

template <typename T, size_t N, size_t M>
bool m::operator==(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) {

    for (size_t x = 0; x < N; x++) {

        for (size_t y = 0; y < M; y++) {

            if (!util::isEqual(lhs.get(x, y), rhs.get(x, y))) return false;
        }
    }

    return true;
}

template <typename T, size_t N, size_t M>
bool m::operator!=(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) {

    return !(lhs == rhs);
}

// TODO: Either make this less crappy or remove it
template <typename T, size_t N, size_t M>
std::ostream &m::operator<<(std::ostream &lhs, const m::tmat<T, N, M> &rhs) {

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

template <typename T>
m::tmat<T, 4, 4> m::mat::scale(const m::tvec<T, 3> &factors) {

    return tmat<T, 4, 4>(factors.get(0), 0,              0,              0,
                         0,              factors.get(1), 0,              0,
                         0,              0,              factors.get(2), 0,
                         0,              0,              0,              1);
}

template <typename T>
m::tmat<T, 4, 4> m::mat::scale(T factor) {

    return scale(tvec<T, 3>(factor, factor, factor));
}

template <typename T>
m::tmat<T, 4, 4> m::mat::translation(const m::tvec<T, 3> &offset) {

    return tmat<T, 4, 4>(1, 0, 0, offset.x(),
                         0, 1, 0, offset.y(),
                         0, 0, 1, offset.z(),
                         0, 0, 0, 1);
}

template <typename T>
m::tmat<T, 4, 4> m::mat::rotation(const m::tquat<T> &rep) {

    tvec<T, 3> rotatedX = rep.rotate(m::X_AXIS<T>);
    tvec<T, 3> rotatedY = rep.rotate(m::Y_AXIS<T>);
    tvec<T, 3> rotatedZ = rep.rotate(m::Z_AXIS<T>);

    return tmat<T, 4, 4>(rotatedX.x(), rotatedY.x(), rotatedZ.x(), 0,
                         rotatedX.y(), rotatedY.y(), rotatedZ.y(), 0,
                         rotatedX.z(), rotatedY.z(), rotatedZ.z(), 0,
                         0,            0,            0,            1);
}

template <typename T>
m::tmat<T, 4, 4> m::mat::rotation(T angle, const m::tvec<T, 3> &axis) {

    return rotation(tquat<T>::rotation(angle, axis));
}

#endif
