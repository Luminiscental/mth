#pragma once

/// @file mth/vec.hpp Header containing definitions for vector expressions.

#ifndef DEFAULT_MEMOIZE
#define DEFAULT_MEMOIZE false
#else
#undef DEFAULT_MEMOIZE
#define DEFAULT_MEMOIZE true
#endif

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace mth
{
    template <typename Derived, bool memoize = DEFAULT_MEMOIZE>
    class tvec_expr;

    template <typename Vec>
    class tvec_info;

    /**
     * @brief Alias for the element type of a given vector expression type.
     *
     * @tparam Vec The vector expression type.
     */
    template <typename Vec>
    using tvec_elem_t = typename tvec_info<Vec>::Elem;

    /**
     * @brief The dimension of a given vector expression type.
     *
     * @tparam Vec The vector expression type.
     */
    template <typename Vec>
    constexpr auto tvec_size_v = tvec_info<Vec>::SIZE;

    template <typename... Ts>
    using enable_if_vecs_t = std::enable_if_t<(
        (std::is_base_of_v<
             tvec_expr<Ts, true>,
             Ts> || std::is_base_of_v<tvec_expr<Ts, false>, Ts>) &&...)>;

    template <typename T>
    using enable_if_vec_t = enable_if_vecs_t<T>;

    /**
     * @brief Vector expression for a concrete / evaluated vector value.
     *
     * This implements the typical iterable protocol to allow range-for loops,
     * whereas the generic `mth::tvec_expr` doesn't; only concrete / evaluated
     * vectors can be iterated. The iterators are forwarded directly from
     * `std::array` for the components.
     *
     * Various aliases for this type are provided with the typical type prefixes
     * and default types: `ivec1` for `tvec<int, 1>`, `dvec3` for `tvec<double,
     * 3>`, `uvec2` for `tvec<unsigned int, 2>`, `fvec4` for `tvec<float, 4>`.
     * Leaving out the type prefix defaults to float, so `vec3` is `fvec3` and
     * so on.
     *
     * This class also implements the compound assignment operators that
     * `tvec_expr` does not since operations on `tvec_expr` change the type and
     * so cannot be used in compound assignment.
     *
     * @tparam T The element type of the vector.
     * @tparam N The dimension of the vector.
     */
    template <typename T, size_t N>
    class tvec : public tvec_expr<tvec<T, N>, false>
    {
    private:
        std::array<T, N> _components;

        template <typename Vec, size_t... Ns>
        constexpr tvec(std::index_sequence<Ns...>, Vec expr)
            : _components{expr[Ns]...}
        {
        }

    public:
        /**
         * @brief Default constructor calls explicit default intialization for
         * all elements.
         *
         * Enabled when the element type `T` is default constructible.
         */
        template <
            typename = std::enable_if_t<std::is_default_constructible_v<T>>>
        constexpr tvec() : _components{}
        {
        }

        /**
         * @brief Convert a vector expression to a concrete vector value.
         *
         * Enabled when `Vec` is derived from `mth::tvec_expr` appropriately.
         *
         * @tparam Vec The vector expression type.
         * @param expr The vector expression value.
         */
        template <
            typename Vec,
            typename = enable_if_vec_t<Vec>,
            typename = std::enable_if_t<
                std::is_same_v<tvec_elem_t<Vec>, T> && tvec_size_v<Vec> == N>>
        constexpr tvec(Vec expr)
            : tvec{std::make_index_sequence<tvec_size_v<Vec>>{},
                   std::move(expr)}
        {
        }

        /**
         * @brief Construct a concrete vector value from a variadic list of N
         * elements of type T.
         *
         * @param values The components of the vector to construct.
         */
        template <
            typename... Ts,
            typename = std::enable_if_t<
                sizeof...(Ts) == N && (std::is_same_v<Ts, T> && ...)>>
        constexpr tvec(Ts... values) : _components{std::move(values)...}
        {
        }

        // Compound assignment operators:

        /**
         * @brief Compound assignment operator for addition, adds the given
         * vector expression to this vector.
         *
         * @tparam Vec The vector expression type to add.
         * @param other The vector expression value to add.
         */
        template <
            typename Vec,
            typename = enable_if_vec_t<Vec>,
            typename = std::enable_if_t<
                std::is_same_v<tvec_elem_t<Vec>, T> && tvec_size_v<Vec> == N>>
        constexpr tvec<T, N>& operator+=(Vec other);

        /**
         * @brief Compound assignment operator for subtraction, subtracts the
         * given vector expression from this vector.
         *
         * @tparam Vec The vector expression type to subtract.
         * @param other The vector expression value to subtract.
         */
        template <
            typename Vec,
            typename = enable_if_vec_t<Vec>,
            typename = std::enable_if_t<
                std::is_same_v<tvec_elem_t<Vec>, T> && tvec_size_v<Vec> == N>>
        constexpr tvec<T, N>& operator-=(Vec other);

        /**
         * @brief Compound assignment operator for scalar multiplication, scales
         * this vector by the given scalar.
         *
         * @param scalar The value to scale by.
         */
        constexpr tvec<T, N>& operator*=(T scalar);

        /**
         * @brief Compound assignment operator for scalar division, divides
         * this vector by the given scalar.
         *
         * @param scalar The value to divide by.
         */
        constexpr tvec<T, N>& operator/=(T scalar);

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The value of the component at the passed index.
         */
        constexpr T operator[](size_t index) const
        {
            return _components[index];
        }

        // iterators:

        constexpr auto begin()
        {
            return _components.begin();
        }

        constexpr auto begin() const
        {
            return _components.begin();
        }

        constexpr auto cbegin() const
        {
            return _components.cbegin();
        }

        constexpr auto rbegin()
        {
            return _components.rbegin();
        }

        constexpr auto crbegin() const
        {
            return _components.crbegin();
        }

        constexpr auto end()
        {
            return _components.end();
        }

        constexpr auto end() const
        {
            return _components.end();
        }

        constexpr auto cend() const
        {
            return _components.cend();
        }

        constexpr auto rend()
        {
            return _components.rend();
        }

        constexpr auto crend() const
        {
            return _components.crend();
        }
    };

    template <typename Vec, typename = enable_if_vec_t<Vec>>
    tvec(Vec expr)->tvec<tvec_elem_t<Vec>, tvec_size_v<Vec>>;

    template <typename T, size_t N>
    class tvec_info<tvec<T, N>>
    {
    public:
        using Elem                 = T;
        static constexpr auto SIZE = N;
    };

    /**
     * @brief Vector expression for a sum of two vectors.
     *
     * @tparam Lhs The vector expression type of the left side operand.
     * @tparam Rhs The vector expression type of the right side operand.
     */
    template <typename Lhs, typename Rhs>
    class tvec_sum : public tvec_expr<tvec_sum<Lhs, Rhs>>
    {
    private:
        Lhs _lhs;
        Rhs _rhs;

    public:
        /**
         * @brief Construct the sum expression from its operands.
         *
         * @param lhs The vector expression value of the left side operand.
         * @param rhs The vector expression value of the right side operand.
         */
        explicit constexpr tvec_sum(Lhs lhs, Rhs rhs)
            : _lhs{std::move(lhs)}, _rhs{std::move(rhs)}
        {
        }

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The sum of the components of the left and right side
         * operands at this index.
         */
        constexpr tvec_elem_t<Lhs> operator[](size_t index) const
        {
            return _lhs[index] + _rhs[index];
        }
    };

    template <typename Lhs, typename Rhs>
    class tvec_info<tvec_sum<Lhs, Rhs>>
    {
    public:
        using Elem                 = tvec_elem_t<Lhs>;
        static constexpr auto SIZE = tvec_size_v<Lhs>;
    };

    /**
     * @brief Vector expression for the difference of two vectors.
     *
     * @tparam Lhs The vector expression type of the left side operand.
     * @tparam Rhs The vector expression type of the right side operand.
     */
    template <typename Lhs, typename Rhs>
    class tvec_diff : public tvec_expr<tvec_diff<Lhs, Rhs>>
    {
    private:
        Lhs _lhs;
        Rhs _rhs;

    public:
        /**
         * @brief Construct the difference expression from its operands.
         *
         * @param lhs The vector expression value of the left side operand.
         * @param rhs The vector expression value of the right side operand.
         */
        explicit constexpr tvec_diff(Lhs lhs, Rhs rhs)
            : _lhs{std::move(lhs)}, _rhs{std::move(rhs)}
        {
        }

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The difference of the components of the left and right side
         * operands at this index.
         */
        constexpr tvec_elem_t<Lhs> operator[](size_t index) const
        {
            return _lhs[index] - _rhs[index];
        }
    };

    template <typename Lhs, typename Rhs>
    class tvec_info<tvec_diff<Lhs, Rhs>>
    {
    public:
        using Elem                 = tvec_elem_t<Lhs>;
        static constexpr auto SIZE = tvec_size_v<Lhs>;
    };

    /**
     * @brief Vector expression for the negation of a vector.
     *
     * @tparam Vec The vector expression type of the operand.
     */
    template <typename Vec>
    class tvec_neg : public tvec_expr<tvec_neg<Vec>>
    {
    private:
        Vec _target;

    public:
        /**
         * @brief Construct the negation expression from its operand value.
         *
         * @param vec The vector expression value of the operand.
         */
        explicit constexpr tvec_neg(Vec vec) : _target{std::move(vec)} {}

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The negation of the component of the operand at this index.
         */
        constexpr tvec_elem_t<Vec> operator[](size_t index) const
        {
            return -_target[index];
        }
    };

    template <typename Vec>
    class tvec_info<tvec_neg<Vec>>
    {
    public:
        using Elem                 = tvec_elem_t<Vec>;
        static constexpr auto SIZE = tvec_size_v<Vec>;
    };

    /**
     * @brief Vector expression for scalar-vector multiplication.
     *
     * @tparam T The scalar type.
     * @tparam Vec The vector expression type of the vector operand.
     */
    template <typename T, typename Vec>
    class tvec_scale : public tvec_expr<tvec_scale<T, Vec>>
    {
    private:
        Vec _vector;
        T _scalar;

    public:
        /**
         * @brief Construct the product expression from its operands.
         *
         * @param scalar The scalar value.
         * @param vector The vector expression value.
         */
        explicit constexpr tvec_scale(T scalar, Vec vector)
            : _vector{std::move(vector)}, _scalar{std::move(scalar)}
        {
        }

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The product of the scalar operand the component of the vector
         * operand at this index.
         */
        constexpr T operator[](size_t index) const
        {
            return _scalar * _vector[index];
        }
    };

    template <typename T, typename Vec>
    class tvec_info<tvec_scale<T, Vec>>
    {
    public:
        using Elem                 = tvec_elem_t<Vec>;
        static constexpr auto SIZE = tvec_size_v<Vec>;
    };

    /**
     * @brief Vector expression for dividing a vector by a scalar.
     *
     * @tparam T The scalar type.
     * @tparam Vec The vector expression type of the vector operand.
     */
    template <typename T, typename Vec>
    class tvec_reduce : public tvec_expr<tvec_reduce<T, Vec>>
    {
    private:
        Vec _vector;
        T _scalar;

    public:
        /**
         * @brief Construct the division expression from its operands.
         *
         * @param scalar The scalar value.
         * @param vector The vector expression value.
         */
        explicit constexpr tvec_reduce(T scalar, Vec vector)
            : _vector{std::move(vector)}, _scalar{std::move(scalar)}
        {
        }

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The component of the vector operand at this index divided by
         * the scalar operand.
         */
        constexpr T operator[](size_t index) const
        {
            return _vector[index] / _scalar;
        }
    };

    template <typename T, typename Vec>
    class tvec_info<tvec_reduce<T, Vec>>
    {
    public:
        using Elem                 = tvec_elem_t<Vec>;
        static constexpr auto SIZE = tvec_size_v<Vec>;
    };

    /**
     * @brief Vector expression for mapping a function over the components of a
     * list of vectors.
     *
     * This is the functional concept of mapping; a function f mapped over
     * vectors `a = {a0, a1, ...}, b = {b0, b1, ...}, c = {c0, c1, ...}` is
     * essentially the vector `{f(a0, b0, c0), f(a1, b1, c1), ...}`.
     *
     * The pack of vectors to map over must be non-empty.
     *
     * @tparam Func The type of the function to map.
     * @tparam Vecs Parameter pack of vector expression types to be mapped over.
     */
    template <typename Func, typename... Vecs>
    class tvec_map : public tvec_expr<tvec_map<Func, Vecs...>>
    {
    private:
        std::tuple<Vecs...> _vecs;
        Func _functor;

        template <size_t... Ns>
        constexpr auto
        indexHelper(std::index_sequence<Ns...>, size_t index) const
        {
            return _functor((std::get<Ns>(_vecs)[index])...);
        }

    public:
        /**
         * @brief Construct the vector expression from the function and vector
         * expression values to map over.
         *
         * @param functor The function value to map.
         * @param vecs The pack of vectors to map over.
         */
        explicit constexpr tvec_map(Func functor, Vecs... vecs)
            : _vecs{std::move(vecs)...}, _functor{std::move(functor)}
        {
        }

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The function applied to the components of the vectors at this
         * index.
         */
        constexpr auto operator[](size_t index) const
        {
            return indexHelper(
                std::make_index_sequence<sizeof...(Vecs)>{}, index);
        }
    };

    // note: specialized to valid parameter packs
    template <typename Func, typename Vec, typename... Vecs>
    class tvec_info<tvec_map<Func, Vec, Vecs...>>
    {
    public:
        using Elem =
            std::result_of_t<Func(tvec_elem_t<Vec>, tvec_elem_t<Vecs>...)>;
        static constexpr auto SIZE = tvec_size_v<Vec>;
    };

    /**
     * @brief Vector expression for casting to a new element type.
     *
     * @tparam T The element type to cast to.
     * @tparam Vec The vector expression type of the vector being cast.
     */
    template <typename T, typename Vec>
    class tvec_cast : public tvec_expr<tvec_cast<T, Vec>>
    {
    private:
        Vec _target;

    public:
        /**
         * @brief Construct the vector expression from the cast target
         * expression.
         *
         * @param target The expression to be cast.
         */
        explicit constexpr tvec_cast(Vec target) : _target{std::move(target)} {}

        /**
         * @brief The `operator[]` to define this as a vector expression.
         *
         * @param index The index of the component to access.
         * @return The component of the target expression at this index
         * statically cast to type `T`.
         */
        constexpr auto operator[](size_t index) const
        {
            return static_cast<T>(_target[index]);
        }
    };

    template <typename T, typename Vec>
    class tvec_info<tvec_cast<T, Vec>>
    {
    public:
        using Elem                 = T;
        static constexpr auto SIZE = tvec_size_v<Vec>;
    };

    /**
     * @brief A namespace to hold "static" functions relating to
     * `mth::tvec_expr`.
     */
    namespace vec
    {
        template <size_t... Ns, typename Lhs, typename Rhs>
        constexpr auto dotHelper(std::index_sequence<Ns...>, Lhs lhs, Rhs rhs);

        /**
         * @brief Calculate the dot product of two vector expressions.
         *
         */
        template <
            typename Lhs,
            typename Rhs,
            typename = enable_if_vecs_t<Lhs, Rhs>,
            typename = std::enable_if_t<tvec_size_v<Lhs> == tvec_size_v<Rhs>>>
        constexpr auto dot(Lhs lhs, Rhs rhs);

        /**
         * @brief Map a function over a list of vector expressions.
         *
         * This is the functional concept of mapping; a function f mapped over
         * vectors `a = {a0, a1, ...}, b = {b0, b1, ...}, c = {c0, c1, ...}` is
         * essentially the vector `{f(a0, b0, c0), f(a1, b1, c1), ...}`.
         *
         * The pack of vectors to map over must be non-empty.
         *
         * @tparam Func The type of the function to map.
         * @tparam Vecs Parameter pack of vector expression types to be mapped
         * over.
         * @param f The function value to map.
         * @param vecs The vectors to map over.
         * @return A vector expression for the result of mapping `f` over
         * `vecs`.
         */
        template <
            typename Func,
            typename... Vecs,
            typename = enable_if_vecs_t<Vecs...>>
        constexpr auto map(Func f, Vecs... vecs);
    }
    /**
     * @brief The base vector expression class, used for the CRTP pattern.
     *
     * This base class contains almost all of the operations on vectors, such as
     * detailed component access and arithmetic operator overloads, equality
     * operator overloads. However, it does not provide iterator member
     * functions or compound assignment operators, those are only accessible
     * from a concrete `mth::tvec` value.
     *
     * @tparam Derived The type of the class deriving from this base for CRTP.
     * @tparam memoize Bool for whether the base class should handle memoizing
     * `operator[]`.
     */
    template <typename Derived, bool memoize>
    class tvec_expr
    {
    public:
        using Elem                 = tvec_elem_t<Derived>;
        static constexpr auto SIZE = tvec_size_v<Derived>;

    protected:
        // CRTP safeguard
        tvec_expr()  = default;
        ~tvec_expr() = default;

        // memoized elements
        mutable std::unordered_map<size_t, Elem> _memos;

        template <size_t... Ns>
        constexpr Elem normHelper(std::index_sequence<Ns...>) const
        {
            return ((component<Ns>() * component<Ns>()) + ...);
        }

    public:
        // Casts:

        /**
         * @brief Operator overload wrapping the `tvec_cast` expression.
         *
         * @return A `tvec_cast` expression targeting this.
         */
        template <typename T>
        constexpr operator tvec_cast<T, Derived>() const
        {
            Derived copy = static_cast<Derived const&>(*this);
            return tvec_cast<T, Derived>{std::move(copy)};
        }

        /**
         * @brief Operator overload for casting to a concrete vector of
         * different element type.
         *
         * @tparam T The element type to cast to.
         * @return A concrete evaluation of the respective `tvec_cast`
         * expression.
         */
        template <
            typename T,
            typename = std::enable_if_t<!std::is_same_v<T, Elem>>>
        constexpr operator tvec<T, SIZE>() const
        {
            return tvec{static_cast<tvec_cast<T, Derived>>(*this)};
        }

        // Component accessors:

        /**
         * @brief The base element accessor, delegating to derived classes by
         * CRTP. This should not implement bounds checking.
         *
         * If `memoize` is true this handles memoizing.
         *
         * @param index The index of the component to access.
         * @return The value of the component at this index.
         */
        constexpr Elem operator[](size_t index) const
        {
            if constexpr (memoize)
            {
                if (_memos.count(index) == 0)
                {
                    auto value = static_cast<Derived const&>(*this)[index];
                    _memos.insert_or_assign(index, value);
                    return value;
                }
                else
                {
                    return _memos.find(index)->second;
                }
            }
            else
            {
                return static_cast<Derived const&>(*this)[index];
            }
        }

        /**
         * @brief A dynamic bounds checking wrapper for `operator[]`.
         *
         * If `index` is too large `std::out_of_range` is thrown.
         *
         * @param index The index of the component to access.
         * @return The value of the component at this index.
         */
        Elem get(size_t index) const
        {
            if (index >= SIZE)
            {
                throw std::out_of_range{
                    "index " + std::to_string(index)
                    + " out of range for vector of dimension "
                    + std::to_string(SIZE)};
            }

            return (*this)[index];
        }

        /**
         * @brief A static component accessor wrapping `operator[]`.
         *
         * By using a template parameter for the index the bounds check can use
         * `static_assert` instead of exceptions.
         *
         * @tparam N The index of the component to access.
         * @return The value of the component at this index.
         */
        template <size_t N>
        constexpr Elem component() const
        {
            static_assert(N < SIZE, "component accessed out of range");
            return (*this)[N];
        }

        // aliases:

        /// @brief alias for component<0>.
        constexpr Elem x() const
        {
            return component<0>();
        }

        /// @brief alias for component<1>.
        constexpr Elem y() const
        {
            return component<1>();
        }

        /// @brief alias for component<2>.
        constexpr Elem z() const
        {
            return component<2>();
        }

        /// @brief alias for component<3>.
        constexpr Elem w() const
        {
            return component<3>();
        }

        /// @brief alias for component<0>.
        constexpr Elem r() const
        {
            return component<0>();
        }

        /// @brief alias for component<1>.
        constexpr Elem g() const
        {
            return component<1>();
        }

        /// @brief alias for component<2>.
        constexpr Elem b() const
        {
            return component<2>();
        }

        /// @brief alias for component<3>.
        constexpr Elem a() const
        {
            return component<3>();
        }

        // transforms / calculations:

        /**
         * @brief Calculate the square magnitude of this vector. The result
         * stays in the type of the elements.
         *
         * @return The square magnitude of the vector expression.
         */
        constexpr Elem norm() const
        {
            return normHelper(std::make_index_sequence<SIZE>{});
        }

        /**
         * @brief Calculate the magnitude (or length) of this vector.
         *
         * The calculation is done as a double for generality, directly
         * equivalent to `std::sqrt(static_cast<double>(norm()))`.
         *
         * @return The magnitude (or length) of this vector.
         */
        constexpr double magn() const
        {
            return std::sqrt(static_cast<double>(norm()));
        }

        /**
         * @brief Calculate the dot product with another vector expression.
         * Delegates to `mth::vec::dot`.
         *
         * @tparam Vec The vector expression type of the vector to dot with.
         * @param other The vector to dot with.
         * @return The dot product of this with other.
         */
        template <typename Vec>
        constexpr auto dot(Vec other)
        {
            Derived copy = static_cast<Derived const&>(*this);
            return vec::dot(std::move(copy), std::move(other));
        }

        /**
         * @brief Member function specializing `mth::vec::map` to the single
         * parameter case.
         *
         * Maps a function over the components of this vector, in the functional
         * sense. So for a vector `v = {v0, v1, ...}`, `v.map(f) = {f(v0),
         * f(v1), ...}` essentially.
         *
         * @tparam Func The type of the function to map.
         * @param f The function value to map.
         * @return A new vector whose components are the application of `f` to
         * this vector's components.
         */
        template <typename Func>
        constexpr auto map(Func f) const
        {
            Derived copy = static_cast<Derived const&>(*this);
            return tvec_map{std::move(f), std::move(copy)};
        }
    };

    namespace vec
    {
        template <size_t... Ns, typename Lhs, typename Rhs>
        constexpr auto dotHelper(std::index_sequence<Ns...>, Lhs lhs, Rhs rhs)
        {
            return ((lhs.get(Ns) * rhs.get(Ns)) + ...);
        }

        template <typename Lhs, typename Rhs, typename, typename>
        constexpr auto dot(Lhs lhs, Rhs rhs)
        {
            return dotHelper(
                std::make_index_sequence<tvec_size_v<Lhs>>{},
                std::move(lhs),
                std::move(rhs));
        }

        template <typename Func, typename... Vecs, typename>
        constexpr auto map(Func f, Vecs... vecs)
        {
            return tvec_map{std::move(f), std::move(vecs)...};
        }
    }

    template <typename Derived, bool memoize>
    class tvec_info<tvec_expr<Derived, memoize>>
    {
    public:
        using Elem                 = tvec_elem_t<Derived>;
        static constexpr auto SIZE = tvec_size_v<Derived>;
    };

    enum class tvec_tag
    {
        Dummy
    };

    // operator overloads:

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vecs_t<Lhs, Rhs>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator+(Lhs lhs, Rhs rhs)
    {
        return tvec_sum{std::move(lhs), std::move(rhs)};
    }

    template <
        typename Vec,
        typename = enable_if_vec_t<Vec>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator-(Vec vec)
    {
        return tvec_neg{std::move(vec)};
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vecs_t<Lhs, Rhs>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator-(Lhs lhs, Rhs rhs)
    {
        return tvec_diff{std::move(lhs), std::move(rhs)};
    }

    template <
        typename T,
        typename Vec,
        typename = enable_if_vec_t<Vec>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator*(T lhs, Vec rhs)
    {
        return tvec_scale{std::move(lhs), std::move(rhs)};
    }

    template <
        typename T,
        typename Vec,
        typename = enable_if_vec_t<Vec>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator*(Vec lhs, T rhs)
    {
        return rhs * lhs;
    }

    template <
        typename T,
        typename Vec,
        typename = enable_if_vec_t<Vec>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator/(Vec lhs, T rhs)
    {
        return tvec_reduce{std::move(rhs), std::move(lhs)};
    }

    template <typename T, size_t N, size_t... Ns>
    constexpr auto
    equalityHelper(std::index_sequence<Ns...>, tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return ((lhs.get(Ns) == rhs.get(Ns)) && ...);
    }

    template <typename T, size_t N>
    constexpr auto operator==(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return equalityHelper(
            std::make_index_sequence<N>{}, std::move(lhs), std::move(rhs));
    }

    template <typename T, size_t N>
    constexpr auto operator!=(tvec<T, N> lhs, tvec<T, N> rhs)
    {
        return !(lhs == rhs);
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vecs_t<Lhs, Rhs>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator==(Lhs lhs, Rhs rhs)
    {
        return tvec{std::move(lhs)} == tvec{std::move(rhs)};
    }

    template <
        typename Lhs,
        typename Rhs,
        typename = enable_if_vecs_t<Lhs, Rhs>,
        tvec_tag = tvec_tag::Dummy>
    constexpr auto operator!=(Lhs lhs, Rhs rhs)
    {
        return !(lhs == rhs);
    }

    template <
        typename Vec,
        typename = enable_if_vec_t<Vec>,
        tvec_tag = tvec_tag::Dummy>
    std::ostream& operator<<(std::ostream& ostream, Vec vec)
    {
        ostream << std::string{"("} << vec[0];

        for (size_t i = 1; i < tvec_size_v<Vec>; i++)
        {
            ostream << std::string{", "} << vec[i];
        }

        ostream << std::string{")"};

        return ostream;
    }

    // aliases:

    using fvec2 = tvec<float, 2>;
    using fvec3 = tvec<float, 3>;
    using fvec4 = tvec<float, 4>;

    using dvec2 = tvec<double, 2>;
    using dvec3 = tvec<double, 3>;
    using dvec4 = tvec<double, 4>;

    using ivec2 = tvec<int, 2>;
    using ivec3 = tvec<int, 3>;
    using ivec4 = tvec<int, 4>;

    using uvec2 = tvec<unsigned int, 2>;
    using uvec3 = tvec<unsigned int, 3>;
    using uvec4 = tvec<unsigned int, 4>;

    using vec2 = fvec2;
    using vec3 = fvec3;
    using vec4 = fvec4;

    // Compound assignment operator implementations:

    template <typename T, size_t N>
    template <typename Vec, typename, typename>
    constexpr tvec<T, N>& tvec<T, N>::operator+=(Vec other)
    {
        *this = *this + other;
        return *this;
    }

    template <typename T, size_t N>
    template <typename Vec, typename, typename>
    constexpr tvec<T, N>& tvec<T, N>::operator-=(Vec other)
    {
        *this = *this - other;
        return *this;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N>& tvec<T, N>::operator*=(T scalar)
    {
        *this = *this * scalar;
        return *this;
    }

    template <typename T, size_t N>
    constexpr tvec<T, N>& tvec<T, N>::operator/=(T scalar)
    {
        *this = *this / scalar;
        return *this;
    }
}
