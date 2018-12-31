#ifndef __m_mat_h__
#define __m_mat_h__

#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>
#include <type_traits>
#include <algorithm>
#include <numeric>

#ifndef __m_vec_h__
#include <m/vec.h>
#endif

#ifndef __m_constants_h__
#include <m/constants.h>
#endif

#ifndef __m_quat_h__
#include <m/quat.h>
#endif

#define __m_mat_content_toggle__ // NOTE: Safeguard so mat_content.h doesn't get included anywhere else

namespace m {

    template <typename T, size_t N, size_t M>
    class tmat;

    // NOTE: Forward declaration of tmat_aug so that it (and its members) can be used in tmat::inverse
    // TODO: Augment non-square matrices?

    template <typename T, size_t N, typename A>
    class tmat_aug {

    private:

        std::array<A, N> aux;
        tmat<T, N, N> matrix;

        void eliminateFromBelow(size_t x, size_t y);
        void eliminateFromRight(size_t x, size_t y);

    public:

        tmat_aug() noexcept;
        tmat_aug(const tmat<T, N, N> &matrix, const std::array<A, N> &aux) noexcept;

        auto coefficients() const; 
        auto auxilary() const;

        // TODO: degenerate solution cases
        auto solve() const;

        // NOTE: returns N on empty rows
        auto leadingIndex(size_t row) const;
        // NOTE: returns zero on empty rows
        T leadingValue(size_t row) const;

        bool columnIsZero(size_t x) const;
        bool rowIsZero(size_t y) const;
        bool hasZeroRow() const;
        bool singular() const;

        void swapRows(size_t a, size_t b); 
        void scaleRow(size_t index, T scalar);
        void addRow(size_t targetRow, size_t sourceRow, T scalar = 1); 

        void setRow(size_t index, const tvec<T, N> &val, A auxVal); 

        auto ordered() const;
        auto rowEchelon() const; 
        auto reducedRowEchelon() const; 

        friend auto &operator<<(std::ostream &lhs, const tmat_aug<T, N, A> &rhs) {

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
    };   

    // NOTE: Implement tmat<2, 2>

#define M 2
#define N 2

    template <typename T>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N
#undef M

    // NOTE: Implement tmat<2, M>

#define N 2

    template <typename T, size_t M>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef N

    // NOTE: Implement tmat<N, 2>

#define M 2

    template <typename T, size_t N>
    class tmat<T, N, M> {

#include <m/mat_content.h>

    };

#undef M

    // NOTE: Implement tmat<N, M>

    template <typename T, size_t N, size_t M>
    class tmat {

#include <m/mat_content.h>

    };

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A>::tmat_aug() noexcept {}

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A>::tmat_aug(const tmat<T, N, N> &matrix, const std::array<A, N> &aux) noexcept
        :matrix(matrix), aux(aux) {}

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::coefficients() const {

        return matrix;
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::auxilary() const {

        return aux;
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::solve() const {

        if (singular()) throw std::invalid_argument("m::exception: solve() called on singular system");

        return reducedRowEchelon().auxilary();
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::leadingIndex(size_t row) const {

        for (size_t i = 0; i < N; i++) {

            if (!util::checkZero(matrix.get(i, row))) return i;
        }

        return N;
    }

    template <typename T, size_t N, typename A>
    inline T tmat_aug<T, N, A>::leadingValue(size_t row) const {

        auto index = leadingIndex(row);

        if (index == N) return 0;

        return matrix.get(index, row);
    }

    template <typename T, size_t N, typename A>
    inline bool tmat_aug<T, N, A>::columnIsZero(size_t x) const {

        if (x > N - 1) throw std::out_of_range("m::exception: column index out of bounds");
        
        for (size_t y = 0; y < N; y++) {

            if (!util::checkZero(matrix.get(x, y))) return false;
        }

        return true;
    }

    template <typename T, size_t N, typename A>
    inline bool tmat_aug<T, N, A>::rowIsZero(size_t y) const {

        if (y > N - 1) throw std::out_of_range("m::exception: row index out of bounds");

        for (size_t x = 0; x < N; x++) {

            if (!util::checkZero(matrix.get(x, y))) return false;
        }

        return true;
    }

    template <typename T, size_t N, typename A>
    inline bool tmat_aug<T, N, A>::hasZeroRow() const {

        for (size_t y = 0; y < N; y++) {

            if (rowIsZero(y)) return true;
        }

        return false;
    }

    template <typename T, size_t N, typename A>
    inline bool tmat_aug<T, N, A>::singular() const {

        return rowEchelon().hasZeroRow();
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::swapRows(size_t a, size_t b) {

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
    inline void tmat_aug<T, N, A>::scaleRow(size_t index, T scalar) {

        matrix.setRow(index, scalar * matrix.getRow(index));
        aux[index] = scalar * aux[index];
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::addRow(size_t targetRow, size_t sourceRow, T scalar) {

        auto addValue = scalar * matrix.getRow(sourceRow);
        auto resultValue = matrix.getRow(targetRow) + addValue;

        matrix.setRow(targetRow, resultValue);

        aux[targetRow] += scalar * aux[sourceRow];
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::setRow(size_t index, const tvec<T, N> &val, A auxVal) {

        matrix.setRow(index, val);
        aux[index] = auxVal;
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::eliminateFromBelow(size_t x, size_t y) {

        auto targetValue = -matrix.get(x, y);

        if (util::checkZero(targetValue)) return;

        for (size_t iy = y + 1; iy < N; iy++) {

            if (leadingIndex(iy) < x) continue;

            auto val = matrix.get(x, iy);

            if (!util::checkZero(val)) {

                addRow(y, iy, targetValue / val); 
                return;
            }
        }

        throw std::invalid_argument("m::exception: eliminateFromBelow() called on non-eliminable element");
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::eliminateFromRight(size_t x, size_t y) {

        auto targetValue = -matrix.get(x, y);

        if (util::checkZero(targetValue)) return;

        for (size_t iy = 0; iy < N; iy++) {

            if (iy == y) continue;

            auto value = matrix.get(x, iy);

            if (!util::checkZero(value) && leadingIndex(iy) >= x) {

                addRow(y, iy, targetValue / value);
                return;
            }
        }

        throw std::invalid_argument("m::exception: eliminateFromRight() called on non-eliminable element");
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::ordered() const {

        using std::begin;
        using std::end;

        std::array<size_t, N> rowIndices;
        std::iota(begin(rowIndices), end(rowIndices), 0);

        auto compareLeadingIndex = [=] (auto a, auto b) -> bool { return leadingIndex(a) < leadingIndex(b); };
        std::sort(begin(rowIndices), end(rowIndices), compareLeadingIndex);

        tmat_aug<T, N, A> result{matrix, aux};

        for (size_t i = 0; i < N; i++) {

            result.setRow(i, matrix.getRow(rowIndices[i]), aux[rowIndices[i]]);
        }

        return result;
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::rowEchelon() const {

        auto result = ordered();

        for (size_t x = 0; x < N - 1; x++) { 

            if (result.columnIsZero(x)) continue;

            for (size_t y = x + 1; y < N; y++) {

                if (result.rowIsZero(y)) break;

                if (util::checkZero(result.matrix.get(x, y))) continue;

                result.eliminateFromRight(x, y);
                result = result.ordered();
                y = x;
            }
        }

        return result.ordered(); 
    }

    template <typename T, size_t N, typename A>
    inline auto tmat_aug<T, N, A>::reducedRowEchelon() const {

        auto result = rowEchelon();

        if (result.hasZeroRow()) throw std::invalid_argument("m::exception: reducedRowEchelon() called on singular matrix");

        for (size_t y = 0; y < N; y++) { 

            auto leadingValue = result.leadingValue(y);

            result.scaleRow(y, 1 / leadingValue);

            for (size_t x = result.leadingIndex(y) + 1; x < N; x++) {

                result.eliminateFromBelow(x, y);
            }
        }

        return result;
    }

    // NOTE: Static functions for matrices 
    
    namespace mat {

        template <typename T>
        auto scale(const tvec<T, 3> &factors) {

            return tmat<T, 4, 4>(factors.get(0), 0,              0,              0,
                                 0,              factors.get(1), 0,              0,
                                 0,              0,              factors.get(2), 0,
                                 0,              0,              0,              1);
        }

        template <typename T>
        auto scale(T factor) {

            return scale(tvec<T, 3>(factor, factor, factor));
        }

        template <typename T>
        auto translate(const tvec<T, 3> &offset) {

            return tmat<T, 4, 4>(1, 0, 0, offset.x(),
                                 0, 1, 0, offset.y(),
                                 0, 0, 1, offset.z(),
                                 0, 0, 0, 1);
        }

        template <typename T>
        auto rotate(const tquat<T> &rep) {

            tvec<T, 3> rotatedX = rep.rotate(m::X_AXIS<T>);
            tvec<T, 3> rotatedY = rep.rotate(m::Y_AXIS<T>);
            tvec<T, 3> rotatedZ = rep.rotate(m::Z_AXIS<T>);

            return tmat<T, 4, 4>(rotatedX.x(), rotatedY.x(), rotatedZ.x(), 0,
                                 rotatedX.y(), rotatedY.y(), rotatedZ.y(), 0,
                                 rotatedX.z(), rotatedY.z(), rotatedZ.z(), 0,
                                 0,            0,            0,            1);
        }

        template <typename T>
        auto rotate(T angle, const tvec<T, 3> &axis) {

            return rotate(tquat<T>::rotation(angle, axis));
        }

        // TODO: Ortho and perspective projections
    }

#define CREATE_ALIASES(n) using imat ## n = tmat<int, n, n>; \
                          using lmat ## n = tmat<long, n, n>; \
                          using  mat ## n = tmat<float, n, n>; \
                          using dmat ## n = tmat<double, n, n>; \
                          using cmat ## n = tmat<std::complex<double>, n, n>;

    CREATE_ALIASES(2)
    CREATE_ALIASES(3)
    CREATE_ALIASES(4)
    CREATE_ALIASES(5)
    CREATE_ALIASES(6)
    CREATE_ALIASES(7)
    CREATE_ALIASES(8)
    CREATE_ALIASES(9)

#undef CREATE_ALIASES

}

#undef __m_mat_content_toggle__

#endif
