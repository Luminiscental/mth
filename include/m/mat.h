#ifndef __m_mat_h__
#define __m_mat_h__

#include <complex>
#include <cmath>

#include <iostream>
#include <iomanip>

#include <array>
#include <type_traits>
#include <algorithm>

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

    // TODO: Non-square matrices
    
    template <typename T, size_t N>
    class tmat;

    // NOTE: Forward declaration of tmat_aug so that it (and its members) can be used in tmat::inverse

    template <typename T, size_t N, typename A>
    class tmat_aug {

    private:

        tmat<T, N> matrix;
        std::array<A, N> aux;

        void eliminate(size_t x, size_t y);

    public:

        tmat_aug();
        tmat_aug(const tmat<T, N> &matrix, const std::array<A, N> &aux);

        const tmat<T, N> &coefficients() const; 
        const std::array<A, N> &auxilary() const;

        // TODO: degenerate solution cases
        std::array<A, N> solve() const;

        // NOTE: returns N on empty rows
        size_t leadingIndex(size_t row) const;
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

        tmat_aug<T, N, A> ordered() const;
        tmat_aug<T, N, A> rowEchelon() const; 
        tmat_aug<T, N, A> reducedRowEchelon() const; 

        friend std::ostream &operator<<(std::ostream &stream, const tmat_aug<T, N, A> &augMatrix) {

            stream << std::fixed << std::setprecision(

#ifdef m_PRECISION

            m_PRECISION

#else

            2

#endif

            );

            for (size_t y = 0; y < N; y++) {

                stream << "|\t";

                for (size_t x = 0; x < N - 1; x++) {

                    stream << augMatrix.matrix.get(x, y) << "\t";
                }

                stream << augMatrix.matrix.get(N - 1, y) << "\t|\t" << augMatrix.aux[y] << "\t|";

                if (y < N - 1) stream << std::endl;
            }
            
            return stream;
        }
    };   

    // NOTE: tmat<T, 2> as base case for minor() template recursion

#define __m_mat_basecaseimpl__

#define N 2

    template <typename T>
    class tmat<T, N> {

#include <m/mat_content.h>

    };

#undef N

#undef __m_mat_basecaseimpl__

    // NOTE: tmat<T, N> for N > 2

    template <typename T, size_t N>
    class tmat {

#include <m/mat_content.h>

    };

    // NOTE: Implementation of tmat_aug

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A>::tmat_aug() {}

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A>::tmat_aug(const tmat<T, N> &matrix, const std::array<A, N> &aux) : matrix(matrix), aux(aux) {}

    template <typename T, size_t N, typename A>
    inline const tmat<T, N> &tmat_aug<T, N, A>::coefficients() const {

        return matrix;
    }

    template <typename T, size_t N, typename A>
    inline const std::array<A, N> &tmat_aug<T, N, A>::auxilary() const {

        return aux;
    }

    template <typename T, size_t N, typename A>
    inline std::array<A, N> tmat_aug<T, N, A>::solve() const {

        if (singular()) throw std::invalid_argument("m::exception: solve() called on singular system");

        return reducedRowEchelon().auxilary();
    }

    template <typename T, size_t N, typename A>
    inline size_t tmat_aug<T, N, A>::leadingIndex(size_t row) const {

        for (size_t i = 0; i < N; i++) {

            if (!util::checkZero(matrix.get(i, row))) return i;
        }

        return N;
    }

    template <typename T, size_t N, typename A>
    inline T tmat_aug<T, N, A>::leadingValue(size_t row) const {

        size_t index = leadingIndex(row);

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

        tvec<T, N> rowA = matrix.getRow(a);
        tvec<T, N> rowB = matrix.getRow(b);

        matrix.setRow(a, rowB);
        matrix.setRow(b, rowA);

        A auxA = aux[a];
        A auxB = aux[b];

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

        tvec<T, N> addValue = scalar * matrix.getRow(sourceRow);
        tvec<T, N> resultValue = matrix.getRow(targetRow) + addValue;

        matrix.setRow(targetRow, resultValue);

        aux[targetRow] += scalar * aux[sourceRow];
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::setRow(size_t index, const tvec<T, N> &val, A auxVal) {

        matrix.setRow(index, val);
        aux[index] = auxVal;
    }

    template <typename T, size_t N, typename A>
    inline void tmat_aug<T, N, A>::eliminate(size_t x, size_t y) {

        T targetValue = -matrix.get(x, y);

        if (util::checkZero(targetValue)) return;

        for (size_t iy = 0; iy < N; iy++) {

            if (iy == y) continue;

            T value = matrix.get(x, iy);

            if (!util::checkZero(value) && leadingIndex(iy) >= x) {

                addRow(y, iy, targetValue / value);
                return;
            }
        }

        throw std::invalid_argument("m::exception: eliminate() called on non-eliminable element");
    }

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A> tmat_aug<T, N, A>::ordered() const {

        std::array<size_t, N> rowIndices;

        for (size_t i = 0; i < N; i++) rowIndices[i] = i;

        std::sort(std::begin(rowIndices), std::end(rowIndices), [=] (auto a, auto b) -> bool { return leadingIndex(a) < leadingIndex(b); });

        tmat_aug<T, N, A> result(matrix, aux);

        for (size_t i = 0; i < N; i++) {

            result.setRow(i, matrix.getRow(rowIndices[i]), aux[rowIndices[i]]);
        }

        return result;
    }

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A> tmat_aug<T, N, A>::rowEchelon() const {

        tmat_aug<T, N, A> result = ordered();

        for (size_t x = 0; x < N - 1; x++) { 

            if (result.columnIsZero(x)) continue;

            for (size_t y = x + 1; y < N; y++) {

                if (result.rowIsZero(y)) break;

                if (util::checkZero(result.matrix.get(x, y))) continue;

                result.eliminate(x, y);
                result = result.ordered();
                y = x;
            }
        }

        return result.ordered(); 
    }

    template <typename T, size_t N, typename A>
    inline tmat_aug<T, N, A> tmat_aug<T, N, A>::reducedRowEchelon() const {

        tmat_aug<T, N, A> result = rowEchelon();

        if (result.hasZeroRow()) throw std::invalid_argument("m::exception: reducedRowEchelon() called on singular matrix");

        for (size_t y = 0; y < N; y++) { 

            T leadingValue = result.leadingValue(y);

            result.scaleRow(y, 1 / leadingValue);

            for (size_t x = result.leadingIndex(y) + 1; x < N; x++) {

                T targetVal = -result.matrix.get(x, y);

                for (size_t iy = y + 1; iy < N; iy++) {

                    if (result.leadingIndex(iy) < x) continue;

                    T val = result.matrix.get(x, iy);

                    if (!util::checkZero(val)) {

                        result.addRow(y, iy, targetVal / val); 
                        break;
                    }
                }
            }
        }

        return result;
    }

    // NOTE: transformations as tmat<T, 4> representations (e.g. for 3d rendering)
    
    namespace mat {

        template <typename T>
        tmat<T, 4> scale(const tvec<T, 3> &factors) {

            return tmat<T, 4>(factors.get(0), 0,              0,              0,
                              0,              factors.get(1), 0,              0,
                              0,              0,              factors.get(2), 0,
                              0,              0,              0,              1);
        }

        template <typename T>
        tmat<T, 4> scale(T factor) {

            return scale(tvec<T, 3>(factor, factor, factor));
        }

        template <typename T>
        tmat<T, 4> translate(const tvec<T, 3> &offset) {

            return tmat<T, 4>(1, 0, 0, offset.x(),
                              0, 1, 0, offset.y(),
                              0, 0, 1, offset.z(),
                              0, 0, 0, 1);
        }

        template <typename T>
        tmat<T, 4> rotate(const tquat<T> &rep) {

            tvec<T, 3> rotatedX = rep.rotate(m::X_AXIS<T>);
            tvec<T, 3> rotatedY = rep.rotate(m::Y_AXIS<T>);
            tvec<T, 3> rotatedZ = rep.rotate(m::Z_AXIS<T>);

            return tmat<T, 4>(rotatedX.x(), rotatedY.x(), rotatedZ.x(), 0,
                              rotatedX.y(), rotatedY.y(), rotatedZ.y(), 0,
                              rotatedX.z(), rotatedY.z(), rotatedZ.z(), 0,
                              0,            0,            0,            1);
        }

        template <typename T>
        tmat<T, 4> rotate(T angle, const tvec<T, 3> &axis) {

            return rotate(tquat<T>::rotation(angle, axis));
        }

        // TODO: Ortho and perspective projections
    }

#define TYPEDEF_MAT(n) typedef tmat<int, n> imat ## n; \
                       typedef tmat<long, n> lmat ## n; \
                       typedef tmat<float, n> mat ## n; \
                       typedef tmat<double, n> dmat ## n; \
                       typedef tmat<std::complex<double>, n> cmat ## n;

    TYPEDEF_MAT(2)
    TYPEDEF_MAT(3)
    TYPEDEF_MAT(4)

#undef TYPEDEF_MAT

}

#undef __m_mat_content_toggle__

#endif
