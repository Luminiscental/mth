#ifndef __m_vec_h__
#define __m_vec_h__

#include <iostream>
#include <iomanip>
#include <array>
#include <type_traits>

#include <m/comp.h>

namespace m {

    template <typename T, size_t N>
    class tvec;

    template <typename T, size_t N>
    auto operator+(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    auto operator-(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    auto operator-(const tvec<T, N> &rhs);

    template <typename T, size_t N>
    auto operator*(const T &lhs, const tvec<T, N> &rhs);

    template <typename T, size_t N>
    auto operator*(const tvec<T, N> &lhs, const T &rhs);

    template <typename T, size_t N>
    auto operator/(const tvec<T, N> &lhs, const T &rhs);

    template <typename T, size_t N>
    auto &operator<<(std::ostream &lhs, const tvec<T, N> &rhs);

    // TODO: noexcept on operators?

    template <typename T, size_t N>
    class tvec {

    private:

        std::array<T, N> values;

    public:

        constexpr tvec() noexcept {

            values.fill(0);
        }

        tvec(const std::array<T, N> &values) noexcept;

        template <typename ...Q, typename std::enable_if<sizeof...(Q) == N, int>::type = 0>
        constexpr tvec(Q... args) noexcept
            :values{static_cast<T>(args)...} {}

        operator std::array<T, N>() noexcept;

        auto begin();

        auto begin() const;

        auto end();

        auto end() const;

        auto cbegin() const;

        auto cend() const;

        auto rbegin();

        auto rend();

        auto crbegin() const;

        auto crend() const;

        constexpr size_t size() const noexcept {

            return N;
        }

        T &get(size_t index);

        const T &get(size_t index) const;

#define BINDING(name, value) const T & name () const; \
                                   T & name ();

        BINDING(x, this->get(0))
        BINDING(y, this->get(1))
        BINDING(z, this->get(2))
        BINDING(w, this->get(3))

#undef BINDING

        // TODO: There should be a better way of doing this

        tvec<T, 2> xy() const;

        tvec<T, 2> yz() const;

        tvec<T, 2> zw() const;
        
        tvec<T, 3> xyz() const;

        tvec<T, 3> yzw() const;

        auto magnSqr() const noexcept;

        auto magn() const noexcept;
        
        auto dot(const tvec<T, N> &rhs) const;

        auto unit() const;

        static auto dot(const tvec<T, N> &lhs, const tvec<T, N> &rhs);

        auto &operator+=(const tvec<T, N> &rhs);

        auto &operator-=(const tvec<T, N> &rhs);

        auto &operator*=(const T &rhs);

        auto &operator/=(const T &rhs);
    };

    // NOTE: Specialized static functions 

    namespace vec {

        template <typename T>
        auto cross(const tvec<T, 3> &lhs, const tvec<T, 3> &rhs);
    }

    // TODO: Maybe move to double as default
    
#define CREATE_ALIASES(n) using ivec ## n = tvec<int, n>; \
                          using lvec ## n = tvec<long, n>; \
                          using  vec ## n = tvec<float, n>; \
                          using dvec ## n = tvec<double, n>; \
                          using cvec ## n = tvec<m::comp, n>;

    CREATE_ALIASES(2)
    CREATE_ALIASES(3)
    CREATE_ALIASES(4)
    CREATE_ALIASES(5)
    CREATE_ALIASES(6)
    CREATE_ALIASES(7)
    CREATE_ALIASES(8)
    CREATE_ALIASES(9)

#undef CREATE_ALIASES

    template <typename T>
    constexpr tvec<T, 3> X_AXIS(1, 0, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Y_AXIS(0, 1, 0);
                                
    template <typename T>
    constexpr tvec<T, 3> Z_AXIS(0, 0, 1);
}

#include <m/vec_impl.h>

#endif
