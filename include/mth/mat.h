#ifndef __mth_mat_h__
#define __mth_mat_h__

/* <mth/mat.h> - matrix header
 *      This includes the tmat template class representing an N by M matrix of coefficients of type T. Addition
 *      and subtraction operations are defined along with multiplication and also division for square matrices.
 *      Member functions to find the determinant and various related constructs such as the inverse and adjoint
 *      are included when N = M.
 * 
 *      This header also includes the tmat_aug template clas which describes a system of linear equations with
 *      N variables and N auxilary values with a separate type. It has defined row operations and methods to solve
 *      the system / reduce to echelon form.
 */

#include <iostream>
#include <array>
#include <type_traits>
#include <numeric>
#include <algorithm>

#include <mth/mth.h>
#include <mth/vec.h>
#include <mth/quat.h>
#include <mth/comp.h>

// Safeguard so <mth/mat_content.h> doesn't get included anywhere else
#define __mth_mat_content_toggle__ 

namespace mth {

    template <typename T, size_t N, size_t M>
    class tmat;

    // Stores an NxN matrix of coefficients of type T with auxilary values of type A
    // Represents the equation matrix . r = aux, where r is a vector with scalar type A
    
    // TODO: This should probably be a Normal class rather than tdata

    template <typename T, size_t N, typename A>
    class tmat_aug {

    private:

        mth::tvec<A, N> aux;

        // TODO: Augment non-square matrices?
        tmat<T, N, N> matrix;

        // Use row operations to set the coefficient at x,y to zero using rows below y without affecting coefficients to the left of it
        void eliminateFromBelow(size_t x, size_t y) {

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
        }

        // Use row operations to set the coefficient at x,y to zero without affecting the coefficients to the left of it
        void eliminateFromRight(size_t x, size_t y) {

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
        }

    public:

        tmat_aug() noexcept = default;
        tmat_aug(const tmat<T, N, N> &matrix, const std::array<A, N> &aux) noexcept
            :matrix(matrix), aux(aux) {}

        tmat<T, N, N> coefficients() const {

            return matrix;
        }

        mth::tvec<A, N> auxilary() const {

            return aux;
        }

        // TODO: degenerate solution cases
        mth::tvec<A, N> solve() const {

            return reducedRowEchelon().auxilary();
        }

        // Index of first non-zero coefficient, returns N on empty rows
        size_t leadingIndex(size_t row) const {

            for (size_t i = 0; i < N; i++) {

                if (!util::isZero(matrix.get(i, row))) return i;
            }

            return N;
        }

        // Value of first non-zero coefficient, returns zero on empty rows
        T leadingValue(size_t row) const {

            auto index = leadingIndex(row);

            if (index == N) return 0;

            return matrix.get(index, row);
        }

        bool columnIsZero(size_t x) const {

            for (size_t y = 0; y < N; y++) {

                if (!util::isZero(matrix.get(x, y))) return false;
            }

            return true;
        }

        bool rowIsZero(size_t y) const {

            for (size_t x = 0; x < N; x++) {

                if (!util::isZero(matrix.get(x, y))) return false;
            }

            return true;
        }

        bool hasZeroRow() const {

            for (size_t y = 0; y < N; y++) {

                if (rowIsZero(y)) return true;
            }

            return false;
        }

        bool singular() const {

            return rowEchelon().hasZeroRow();
        }

        // Row operations

        void swapRows(size_t a, size_t b) {

            auto rowA = matrix.getRow(a);
            auto rowB = matrix.getRow(b);

            matrix.setRow(a, rowB);
            matrix.setRow(b, rowA);

            auto auxA = aux[a];
            auto auxB = aux[b];

            aux[a] = auxB;
            aux[b] = auxA;
        } 

        void scaleRow(size_t index, T scalar) {

            matrix.setRow(index, scalar * matrix.getRow(index));
            aux[index] = scalar * aux[index];
        }

        // Adds sourceRow * scalar to targetRow
        void addRow(size_t targetRow, size_t sourceRow, T scalar = 1) {

            auto addValue = scalar * matrix.getRow(sourceRow);
            auto resultValue = matrix.getRow(targetRow) + addValue;

            matrix.setRow(targetRow, resultValue);

            aux[targetRow] += scalar * aux[sourceRow];
        } 

        void setRow(size_t index, const tvec<T, N> &val, A auxVal) {

            matrix.setRow(index, val);
            aux[index] = auxVal;
        }

        // Transformations

        // Re-orders rows from lowest leadingIndex at the top
        tmat_aug<T, N, A> ordered() const {

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

        // Convert to row echelon form (zero below diagonal)
        tmat_aug<T, N, A> rowEchelon() const {

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

        // Convert to row reduced echelon form (identity matrix)
        tmat_aug<T, N, A> reducedRowEchelon() const {

            // Start from echelon form so that half of the values to eliminate are already zero
            auto result = rowEchelon();

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

        // TODO: Arithmetic on augmented matrices?
    };   

    template <typename T, size_t N, typename A>
    std::ostream &operator<<(std::ostream &lhs, const tmat_aug<T, N, A> &rhs) {

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

    // 0-dimension matrices are undefined
    
    template <typename T, size_t N>
    class tmat<T, N, 0>;

    template <typename T, size_t M>
    class tmat<T, 0, M>;

    template <typename T>
    class tmat<T, 0, 0>;

#define M 1
#define N 1

    // Define tmat<T, 1, 1>

#include <mth/mat_content.h>

#undef N
#undef M

#define N 1

    // Define tmat<T, 1, M> where M > 1

#include <mth/mat_content.h>

#undef N

#define M 1

    // Define tmat<T, N, 1> where N > 1

#include <mth/mat_content.h>

#undef M

    // Define tmat<T, N, M> where N > 1 and M > 1

#include <mth/mat_content.h>

    template <typename T, size_t N, size_t M>
    constexpr tmat<T, N, M> operator+(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) noexcept {

        auto result = lhs;

        return result += rhs;
    }

    template <typename T, size_t N, size_t M>
    constexpr tmat<T, N, M> operator-(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) noexcept {

        auto result = lhs;

        return result -= rhs;
    }

    template <typename T, size_t N, size_t M>
    constexpr tmat<T, N, M> operator*(const T &lhs, const tmat<T, N, M> &rhs) noexcept {

        auto result = rhs;

        return result *= lhs;
    }

    template <typename T, size_t N, size_t M>
    constexpr tmat<T, N, M> operator*(const tmat<T, N, M> &lhs, const T &rhs) {

        return rhs * lhs;
    }

    template <typename T, size_t N, size_t M>
    constexpr tmat<T, N, M> operator/(const tmat<T, N, M> &lhs, const T &rhs) noexcept {

        auto result = lhs;

        return result /= rhs;
    }

    template <typename T, size_t N, size_t M, size_t O>
    constexpr tmat<T, O, M> operator*(const tmat<T, N, M> &lhs, const tmat<T, O, N> &rhs) noexcept {

        tmat<T, O, M> result;

        for (size_t x = 0; x < O; x++) {

            for (size_t y = 0; y < M; y++) {

                result.get(x, y) = tvec<T, N>::dot(lhs.getRow(y), rhs.getColumn(x));
            }
        }

        return result;
    }

    template <typename T, size_t N>
    constexpr tmat<T, N, N> operator/(const tmat<T, N, N> &lhs, const tmat<T, N, N> &rhs) noexcept {

        return lhs * rhs.inverse();
    }

    template <typename T, size_t N, size_t M>
    constexpr tvec<T, M> operator*(const tmat<T, N, M> &lhs, const tvec<T, N> &rhs) noexcept {

        tvec<T, M> result;

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                result.get(y) = tvec<T, N>::dot(lhs.getRow(y), rhs);
            }
        }

        return result;
    }

    template <typename T, size_t N, size_t M>
    constexpr bool operator==(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) noexcept {

        for (size_t x = 0; x < N; x++) {

            for (size_t y = 0; y < M; y++) {

                if (!util::isEqual(lhs.get(x, y), rhs.get(x, y))) return false;
            }
        }

        return true;
    }

    template <typename T, size_t N, size_t M>
    constexpr bool operator!=(const tmat<T, N, M> &lhs, const tmat<T, N, M> &rhs) noexcept {

        return !(lhs == rhs);
    }

    template <typename T, size_t N, size_t M>
    std::ostream &operator<<(std::ostream &lhs, const tmat<T, N, M> &rhs) {

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

    // Specialized "static" functions 
    namespace mat {

        // Transformations assumed to act on (xyz, 1) form 4-vectors

        // Transformation scaling each axis by it's dot with factors
        template <typename T>
        constexpr tmat<T, 4, 4> scale(const tvec<T, 3> &factors) noexcept {

            return tmat<T, 4, 4>(factors.get(0), 0,              0,              0,
                                 0,              factors.get(1), 0,              0,
                                 0,              0,              factors.get(2), 0,
                                 0,              0,              0,              1);
        }

        // Transformation scaling each axis by factor
        template <typename T>
        constexpr tmat<T, 4, 4> scale(T factor) noexcept {

            return scale(tvec<T, 3>(factor, factor, factor));
        }

        template <typename T>
        constexpr tmat<T, 4, 4> translation(const tvec<T, 3> &offset) noexcept {

            return tmat<T, 4, 4>(1, 0, 0, offset.x(),
                                 0, 1, 0, offset.y(),
                                 0, 0, 1, offset.z(),
                                 0, 0, 0, 1);
        }

        // Convert a quaternion rotation representation into its matrix representative
        template <typename T>
        constexpr tmat<T, 4, 4> rotation(const tquat<T> &rep) noexcept {

            tvec<T, 3> rotatedX = rep.rotate(mth::X_AXIS<T>);
            tvec<T, 3> rotatedY = rep.rotate(mth::Y_AXIS<T>);
            tvec<T, 3> rotatedZ = rep.rotate(mth::Z_AXIS<T>);

            return tmat<T, 4, 4>(rotatedX.x(), rotatedY.x(), rotatedZ.x(), 0,
                                 rotatedX.y(), rotatedY.y(), rotatedZ.y(), 0,
                                 rotatedX.z(), rotatedY.z(), rotatedZ.z(), 0,
                                 0,            0,            0,            1);
        }

        // Convert euler angle and axis to its matrix representative
        template <typename T>
        constexpr tmat<T, 4, 4> rotation(T angle, const tvec<T, 3> &axis) noexcept {

            return rotation(tquat<T>::rotation(angle, axis));
        }

        template <typename T>
        constexpr tmat<T, 4, 4> orthographic(T left, T right, T bottom, T top, T near, T far) noexcept {

            auto rml = right - left;
            auto tmb = top - bottom;
            auto fmn = far - near;

            auto rpl = right + left;
            auto tpb = top + bottom;
            auto fpn = far + near;

            return tmat<T, 4, 4>(2 / rml, 0,       0,        -rpl / rml,
                                 0,       2 / tmb, 0,        -tpb / tmb,
                                 0,       0,       -2 / fmn, -fpn / fmn,
                                 0,       0,       0,        1);
        }

        template <typename T>
        constexpr tmat<T, 4, 4> perspective(T left, T right, T bottom, T top, T near, T far) noexcept {

            auto rml = right - left;
            auto tmb = top - bottom;
            auto fmn = far - near;

            auto rpl = right + left;
            auto tpb = top + bottom;
            auto fpn = far + near;

            return tmat<T, 4, 4>(2 * near / rml, 0,              rpl / rml,  0,
                                 0,              2 * near / tmb, tpb / tmb,  0,
                                 0,              0,              -fpn / fmn, -2 * far * near / fmn,
                                 0,              0,              -1,         0);
        }
    }

    // Alias types for single-digit dimensions with scalar types int, long, float, double, mth::comp

#define CREATE_SQR_ALIASES(n)       using imat ## n = tmat<int, n, n>; \
                                    using lmat ## n = tmat<long, n, n>; \
                                    using fmat ## n = tmat<float, n, n>; \
                                    using dmat ## n = tmat<double, n, n>; \
                                    using cmat ## n = tmat<comp, n, n>; \
                                    using  mat ## n = dmat ## n;

#define CREATE_RECT_ALIASES(n, m)   using imat ## n ## x ## m = tmat<int, n, m>; \
                                    using lmat ## n ## x ## m = tmat<long, n, m>; \
                                    using fmat ## n ## x ## m = tmat<float, n, m>; \
                                    using dmat ## n ## x ## m = tmat<double, n, m>; \
                                    using cmat ## n ## x ## m = tmat<comp, n, m>; \
                                    using  mat ## n ## x ## m = dmat ## n ## x ## m;

#define CREATE_ALIASES(n)   CREATE_SQR_ALIASES(n) \
                            CREATE_RECT_ALIASES(n, 1) \
                            CREATE_RECT_ALIASES(n, 2) \
                            CREATE_RECT_ALIASES(n, 3) \
                            CREATE_RECT_ALIASES(n, 4) \
                            CREATE_RECT_ALIASES(n, 5) \
                            CREATE_RECT_ALIASES(n, 6) \
                            CREATE_RECT_ALIASES(n, 7) \
                            CREATE_RECT_ALIASES(n, 8) \
                            CREATE_RECT_ALIASES(n, 9)

    CREATE_ALIASES(1)
    CREATE_ALIASES(2)
    CREATE_ALIASES(3)
    CREATE_ALIASES(4)
    CREATE_ALIASES(5)
    CREATE_ALIASES(6)
    CREATE_ALIASES(7)
    CREATE_ALIASES(8)
    CREATE_ALIASES(9)

#undef CREATE_ALIASES
#undef CREATE_SQR_ALIASES

}

namespace std {

    template <typename T, size_t N, size_t M>
    struct hash<mth::tmat<T, N, M>> {

        size_t operator()(const mth::tmat<T, N, M> &x) {

            size_t result = 0;

            for (auto row : x.rows()) {

                result ^= hash<decltype(row)>()(row);
            }

            return result;
        }
    };
}

#undef __mth_mat_content_toggle__

#endif
