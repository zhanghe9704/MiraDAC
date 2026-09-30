/**
 * @file davector.h
 * @brief DAVector<T>: wrapper around a pool slot with all DA arithmetic.
 *
 * @details Template on coefficient type T.  All methods delegate to the
 *   kernels in da::detail (engine.h).  The DAEnv* is set at construction
 *   time from da_current_env(); copy/move ctors inherit env_ from source.
 *
 *   Every binary DA×DA operation runs the check_env() guard first.
 *
 * Stage 5: numerical complete (T = double).
 */

#pragma once

#include "da/env.h"
#include "da/engine.h"

#include <complex>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <stdexcept>

#ifdef DA_WITH_SYMBOLIC
#  include "da/symbolic_ops.h"
#  include <symengine/expression.h>
#  include <symengine/lambda_double.h>
#  include <symengine/symengine_rcp.h>
#  include <symengine/eval_double.h>
#  include <map>
#  include <functional>
#  include <memory>
#endif

namespace da {

// ===========================================================================
// Forward declarations needed for operators
// ===========================================================================
template<class T> struct DAVector;

using NDA = DAVector<double>;

#ifdef DA_WITH_SYMBOLIC
// SDA alias is defined after the template definition below
#endif

// ===========================================================================
// DAVector<T>
// ===========================================================================

template<class T>
struct DAVector {
    // --------------------------------------------------------------------- //
    //  Data members                                                          //
    // --------------------------------------------------------------------- //
    DAEnv*   env_;   ///< Environment this vector belongs to (never null after ctor)
    unsigned slot_;  ///< Index into env_->template pool<T>()

    // --------------------------------------------------------------------- //
    //  Static epsilon (global threshold for iszero / clean)                 //
    // --------------------------------------------------------------------- //
    static double eps;

    // --------------------------------------------------------------------- //
    //  Constructors                                                          //
    // --------------------------------------------------------------------- //

    /// Default: alloc a zero-const slot from the current env.
    DAVector() : env_(&da_current_env()), slot_(env_->template pool<T>().alloc()) {}

    /// Copy: allocate a new slot and copy data; inherit env_ from source.
    DAVector(const DAVector& other)
        : env_(other.env_), slot_(env_->template pool<T>().alloc())
    {
        detail::ad_copy(env_->layout(), env_->template pool<T>(), other.slot_, slot_);
    }

    /// Marks a vector that owns no slot: the source of a move, or a vector
    /// after clear(). Such a vector may only be destroyed or assigned to.
    static constexpr unsigned no_slot = std::numeric_limits<unsigned>::max();

    /// Move: take over the slot; the source is left owning none, so a move
    /// needs no free slot and cannot throw.
    DAVector(DAVector&& other) noexcept
        : env_(other.env_), slot_(other.slot_)
    {
        other.slot_ = no_slot;
    }

    /// Construct from double constant.
    explicit DAVector(double x)
        : env_(&da_current_env()), slot_(env_->template pool<T>().alloc())
    {
        T val{static_cast<T>(x)};
        detail::ad_const(env_->layout(), env_->template pool<T>(), slot_, val);
    }

    /// Construct from int constant.
    explicit DAVector(int x)
        : env_(&da_current_env()), slot_(env_->template pool<T>().alloc())
    {
        T val{static_cast<T>(static_cast<double>(x))};
        detail::ad_const(env_->layout(), env_->template pool<T>(), slot_, val);
    }

    /// Construct from std::vector<double> of coefficients.
    explicit DAVector(std::vector<double>& v)
        : env_(&da_current_env()), slot_(env_->template pool<T>().alloc())
    {
        Pool<T>& pool = env_->template pool<T>();
        unsigned len = static_cast<unsigned>(
            std::min(v.size(), static_cast<std::size_t>(pool.full_len())));
        T* dst = pool.slot(slot_);
        for (unsigned i = 0; i < len; ++i)
            dst[i] = static_cast<T>(v[i]);
        pool.set_len(slot_, len);
    }

    // --------------------------------------------------------------------- //
    //  Destructor — return the slot to the pool                              //
    // --------------------------------------------------------------------- //
    ~DAVector() {
        // da_clear() (and da_destroy_env()) retire an environment that still
        // has live vectors instead of deleting it: the shell stays valid with
        // its pools released, so the poolsize() check below reads 0 and the
        // free is skipped. That makes "clear, then let vectors go out of
        // scope" safe, as it was in the reference library.
        //
        // Still unsafe: `delete &env` on a DAEnv that vectors reference.
        // Use da_destroy_env() instead.
        if (slot_ != no_slot && env_ && env_->template pool<T>().poolsize() > 0)
            env_->template pool<T>().free(slot_);
    }

    // --------------------------------------------------------------------- //
    //  Assignment operators                                                  //
    // --------------------------------------------------------------------- //

    DAVector& operator=(const DAVector& other) {
        if (this == &other) return *this;
        // Keep our own slot but copy data from other (env_ unchanged).
        if (slot_ == no_slot) slot_ = env_->template pool<T>().alloc();
        detail::ad_copy(env_->layout(), env_->template pool<T>(), other.slot_, slot_);
        return *this;
    }

    DAVector& operator=(DAVector&& other) noexcept {
        if (this == &other) return *this;
        // Swap env and slot together (a slot index only means something in
        // its own env's pool); the other's dtor frees what we give it.
        std::swap(env_, other.env_);
        std::swap(slot_, other.slot_);
        if (other.slot_ != no_slot)
            detail::ad_reset(other.env_->layout(), other.env_->template pool<T>(), other.slot_);
        return *this;
    }

    DAVector& operator=(double x) {
        reset_const(static_cast<T>(x));
        return *this;
    }

    DAVector& operator=(int x) {
        reset_const(static_cast<T>(static_cast<double>(x)));
        return *this;
    }

    // --------------------------------------------------------------------- //
    //  Compound assignment                                                   //
    // --------------------------------------------------------------------- //

    DAVector& operator+=(const DAVector& other) {
        check_env(env_, other.env_);
        detail::ad_add(env_->layout(), env_->template pool<T>(), slot_, other.slot_);
        return *this;
    }

    DAVector& operator+=(DAVector&& other) {
        *this = *this + other;
        return *this;
    }

    DAVector& operator+=(double x) {
        detail::ad_add_const(env_->layout(), env_->template pool<T>(), slot_, static_cast<T>(x));
        return *this;
    }

    DAVector& operator+=(int x) {
        detail::ad_add_const(env_->layout(), env_->template pool<T>(), slot_,
                             static_cast<T>(static_cast<double>(x)));
        return *this;
    }

    DAVector& operator-=(const DAVector& other) {
        check_env(env_, other.env_);
        detail::ad_sub(env_->layout(), env_->template pool<T>(), slot_, other.slot_);
        return *this;
    }

    DAVector& operator-=(DAVector&& other) {
        *this = *this - other;
        return *this;
    }

    DAVector& operator-=(double x) {
        detail::ad_add_const(env_->layout(), env_->template pool<T>(), slot_, static_cast<T>(-x));
        return *this;
    }

    DAVector& operator-=(int x) {
        detail::ad_add_const(env_->layout(), env_->template pool<T>(), slot_,
                             static_cast<T>(static_cast<double>(-x)));
        return *this;
    }

    DAVector& operator*=(const DAVector& other) {
        *this = *this * other;
        return *this;
    }

    DAVector& operator*=(DAVector&& other) {
        *this = *this * other;
        return *this;
    }

    DAVector& operator*=(double x) {
        detail::ad_mult_const(env_->layout(), env_->template pool<T>(), slot_, static_cast<T>(x));
        return *this;
    }

    DAVector& operator*=(int x) {
        detail::ad_mult_const(env_->layout(), env_->template pool<T>(), slot_,
                              static_cast<T>(static_cast<double>(x)));
        return *this;
    }

    DAVector& operator/=(const DAVector& other) {
        *this = *this / other;
        return *this;
    }

    DAVector& operator/=(DAVector&& other) {
        *this = *this / other;
        return *this;
    }

    DAVector& operator/=(double x) {
        detail::ad_div_c(env_->layout(), env_->template pool<T>(), slot_, static_cast<T>(x));
        return *this;
    }

    DAVector& operator/=(int x) {
        detail::ad_div_c(env_->layout(), env_->template pool<T>(), slot_,
                         static_cast<T>(static_cast<double>(x)));
        return *this;
    }

    // --------------------------------------------------------------------- //
    //  Output                                                                //
    // --------------------------------------------------------------------- //

    void print() const;

    // --------------------------------------------------------------------- //
    //  Queries                                                               //
    // --------------------------------------------------------------------- //

    /// Return the constant term (index 0).
    T con() const {
        const T* p = env_->template pool<T>().slot(slot_);
        return p[0];
    }

    /// Number of stored terms (current length of the slot).
    unsigned int length() const {
        return env_->template pool<T>().len(slot_);
    }

    /// Number of non-zero coefficients.
    int n_element() const {
        return detail::ad_n_element(env_->layout(), env_->template pool<T>(), slot_);
    }

    // --------------------------------------------------------------------- //
    //  Element access                                                        //
    // --------------------------------------------------------------------- //

    /// Return the i-th element value (0-based) and its exponent vector c.
    void element(unsigned int i, unsigned int* c, double& elem) const {
        unsigned int ii = i + 1;
        detail::ad_elem(env_->layout(), env_->template pool<T>(), slot_, ii, c,
                        reinterpret_cast<T&>(elem));
    }

    void element(unsigned int i, std::vector<unsigned int>& c, double& elem) const {
        unsigned int ii = i + 1;
        if (c.size() != static_cast<std::size_t>(env_->layout().num_vars()))
            c.resize(env_->layout().num_vars());
        T val{};
        detail::ad_elem(env_->layout(), env_->template pool<T>(), slot_, ii, c.data(), val);
        elem = static_cast<double>(val);
    }

    void derivative(unsigned int i, unsigned int* c, double& elem) const {
        element(i, c, elem);
    }

    void derivative(unsigned int i, std::vector<unsigned int>& c, double& elem) const {
        element(i, c, elem);
    }

    /// Return the scalar value of the i-th element (0-based).
    double element(int i) {
        unsigned int ii = static_cast<unsigned int>(i + 1);
        int d = dim();
        std::vector<unsigned int> c(d);
        T val{};
        detail::ad_elem(env_->layout(), env_->template pool<T>(), slot_, ii, c.data(), val);
        return static_cast<double>(val);
    }

    /// Return the scalar value of an element by exponent vector.
    double element(std::vector<int> idx) {
        T val{};
        detail::ad_pek(env_->layout(), env_->template pool<T>(), slot_,
                       idx.data(), idx.size(), val);
        return static_cast<double>(val);
    }

    /// Return partial derivative value by exponent vector.
    double derivative(std::vector<int> idx) {
        double val = element(idx);
        double fac = 1.0;
        for (int o : idx) {
            for (int k = 1; k <= o; ++k) fac *= k;
        }
        return val * fac;
    }

    /// Return the exponent vector of the i-th element (uses Layout order table).
    const std::vector<int>& element_orders(int i) const {
        return env_->layout().orders_ref(i);
    }

    // --------------------------------------------------------------------- //
    //  Norm                                                                  //
    // --------------------------------------------------------------------- //

    double norm() const {
        return detail::ad_norm(env_->layout(), env_->template pool<T>(), slot_);
    }

    double weighted_norm(double w) const {
        return detail::ad_weighted_norm(env_->layout(), env_->template pool<T>(), slot_, w);
    }

    // --------------------------------------------------------------------- //
    //  Modification                                                          //
    // --------------------------------------------------------------------- //

    /// Set the element at exponent vector c to elem.
    void set_element(int* c, double elem) {
        std::size_t n = static_cast<std::size_t>(env_->layout().num_vars());
        T val{static_cast<T>(elem)};
        detail::ad_pok(env_->layout(), env_->template pool<T>(), slot_, c, n, val);
    }

    void set_element(std::vector<int> idx, double elem) {
        T val{static_cast<T>(elem)};
        detail::ad_pok(env_->layout(), env_->template pool<T>(), slot_,
                       idx.data(), idx.size(), val);
    }

    /// Reset all coefficients to zero (length unchanged).
    void reset() {
        detail::ad_reset_vector(env_->layout(), env_->template pool<T>(), slot_);
    }

    /// Reset to a constant; length = 1.
    void reset_const(T x = T{}) {
        detail::ad_const(env_->layout(), env_->template pool<T>(), slot_, x);
    }

    /// Set small coefficients to zero.
    void clean(double threshold) {
        detail::ad_clean(env_->layout(), env_->template pool<T>(), slot_, threshold);
    }

    void clean() {
        detail::ad_clean(env_->layout(), env_->template pool<T>(), slot_, eps);
    }

    /// Return true if all |coeff| < eps.
    bool iszero() const {
        return detail::ad_zero_check(env_->layout(), env_->template pool<T>(), slot_, eps);
    }

    bool iszero(double threshold) const {
        return detail::ad_zero_check(env_->layout(), env_->template pool<T>(), slot_, threshold);
    }

    /// Return slot to pool. The vector then owns no slot: it may only be
    /// destroyed or assigned to.
    void clear() {
        if (slot_ == no_slot) return;
        env_->template pool<T>().free(slot_);
        slot_ = no_slot;
    }

    // --------------------------------------------------------------------- //
    //  Copy to std::vector                                                   //
    // --------------------------------------------------------------------- //

    void to_vector(std::vector<double>& v) const {
        int len = full_length();
        v.resize(static_cast<std::size_t>(len));
        const T* src = env_->template pool<T>().slot(slot_);
        for (int i = 0; i < len; ++i)
            v[static_cast<std::size_t>(i)] = static_cast<double>(src[i]);
    }

    void to_vector(int length, std::vector<double>& v) const {
        int fl = full_length();
        if (length > fl) length = fl;
        v.resize(static_cast<std::size_t>(length));
        const T* src = env_->template pool<T>().slot(slot_);
        for (int i = 0; i < length; ++i)
            v[static_cast<std::size_t>(i)] = static_cast<double>(src[i]);
    }

    // --------------------------------------------------------------------- //
    //  Static environment queries                                            //
    // --------------------------------------------------------------------- //

    static int dim() {
        return static_cast<int>(da_current_env().layout().num_vars());
    }

    static int order() {
        return static_cast<int>(da_current_env().layout().max_order());
    }

    static int full_length() {
        return static_cast<int>(da_current_env().layout().full_len());
    }

// ===========================================================================
// Symbolic API — all members below are guarded by DA_WITH_SYMBOLIC
// They are defined for all T but use if constexpr so they only work for
// T = SymEngine::Expression.  This allows explicit instantiation of
// DAVector<double> without compile errors.
// ===========================================================================
#ifdef DA_WITH_SYMBOLIC

    // Construct from SymEngine::Expression constant.
    // Only meaningful/usable for T = Expression.
    explicit DAVector(SymEngine::Expression s)
        : env_(&da_current_env()),
          slot_([]()->unsigned{
              if constexpr (std::is_same_v<T, SymEngine::Expression>)
                  return da_current_env().template pool<SymEngine::Expression>().alloc();
              else
                  return da_current_env().template pool<double>().alloc();
          }())
    {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            detail::ad_const(env_->layout(),
                             env_->template pool<SymEngine::Expression>(),
                             slot_, s);
        }
    }

    // Construct from std::vector<Expression> of coefficients.
    explicit DAVector(std::vector<SymEngine::Expression>& v)
        : env_(&da_current_env()),
          slot_([]()->unsigned{
              if constexpr (std::is_same_v<T, SymEngine::Expression>)
                  return da_current_env().template pool<SymEngine::Expression>().alloc();
              else
                  return da_current_env().template pool<double>().alloc();
          }())
    {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            Pool<SymEngine::Expression>& pool = env_->template pool<SymEngine::Expression>();
            unsigned len = static_cast<unsigned>(
                std::min(v.size(), static_cast<std::size_t>(pool.full_len())));
            SymEngine::Expression* dst = pool.slot(slot_);
            for (unsigned i = 0; i < len; ++i)
                dst[i] = v[i];
            pool.set_len(slot_, len);
        }
    }

    // Assignment from Expression
    DAVector& operator=(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            detail::ad_const(env_->layout(),
                             env_->template pool<SymEngine::Expression>(),
                             slot_, x);
        }
        return *this;
    }

    // Compound assignments with Expression — only valid for T=Expression
    DAVector& operator+=(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_add_const(env_->layout(),
                                 env_->template pool<SymEngine::Expression>(), slot_, x);
        return *this;
    }
    DAVector& operator-=(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_add_const(env_->layout(),
                                 env_->template pool<SymEngine::Expression>(), slot_,
                                 SymEngine::Expression(-1) * x);
        return *this;
    }
    DAVector& operator*=(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_mult_const(env_->layout(),
                                  env_->template pool<SymEngine::Expression>(), slot_, x);
        return *this;
    }
    DAVector& operator/=(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_div_c(env_->layout(),
                             env_->template pool<SymEngine::Expression>(), slot_, x);
        return *this;
    }

    // element_expr: return Expression value (0-based index)
    SymEngine::Expression element_expr(int i) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            unsigned int ii = static_cast<unsigned int>(i + 1);
            int d = dim();
            std::vector<unsigned int> c(static_cast<std::size_t>(d));
            SymEngine::Expression val{};
            detail::ad_elem(env_->layout(),
                            env_->template pool<SymEngine::Expression>(),
                            slot_, ii, c.data(), val);
            return val;
        } else {
            return SymEngine::Expression(0);
        }
    }

    // element_expr by exponent vector
    SymEngine::Expression element_expr(std::vector<int> idx) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            SymEngine::Expression val{};
            detail::ad_pek(env_->layout(),
                           env_->template pool<SymEngine::Expression>(),
                           slot_, idx.data(), idx.size(), val);
            return val;
        } else {
            return SymEngine::Expression(0);
        }
    }

    // set_element with Expression
    void set_element(int* c, SymEngine::Expression elem) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            std::size_t n = static_cast<std::size_t>(env_->layout().num_vars());
            detail::ad_pok(env_->layout(),
                           env_->template pool<SymEngine::Expression>(), slot_, c, n, elem);
        }
    }

    void set_element(std::vector<int> idx, SymEngine::Expression elem) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_pok(env_->layout(),
                           env_->template pool<SymEngine::Expression>(),
                           slot_, idx.data(), idx.size(), elem);
    }

    void set_element(int i, SymEngine::Expression elem) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            std::vector<int> idx_vec(env_->layout().orders_ref(i).begin(),
                                     env_->layout().orders_ref(i).end());
            set_element(idx_vec, elem);
        }
    }

    // reset_const with Expression — only valid/meaningful for T=Expression.
    // Named reset_const_expr to avoid overload conflict with reset_const(T).
    void reset_const_expr(SymEngine::Expression x) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>)
            detail::ad_const(env_->layout(),
                             env_->template pool<SymEngine::Expression>(),
                             slot_, x);
    }

    // Simplify all symbolic coefficients in-place
    void simplify() {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            unsigned len = env_->template pool<SymEngine::Expression>().len(slot_);
            SymEngine::Expression* p =
                env_->template pool<SymEngine::Expression>().slot(slot_);
            for (unsigned i = 0; i < len; ++i) {
                if (!da::is_zero(p[i]))
                    da::simplified_expr(p[i]);
            }
        }
    }

    // eval: substitute symbol values in-place, converting to numeric SDA
    void eval(SymEngine::map_basic_basic m) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            unsigned len = env_->template pool<SymEngine::Expression>().len(slot_);
            SymEngine::Expression* p =
                env_->template pool<SymEngine::Expression>().slot(slot_);
            for (unsigned i = 0; i < len; ++i) {
                SymEngine::Expression subs_expr = p[i].subs(m);
                p[i] = SymEngine::Expression(
                    SymEngine::eval_double(*subs_expr.get_basic()));
            }
        }
    }

    // eval_funs: produce a vector of LambdaRealDoubleVisitors (one per element)
    void eval_funs(std::vector<SymEngine::RCP<const SymEngine::Basic>> vars,
                   std::vector<SymEngine::LambdaRealDoubleVisitor>& vec) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            unsigned len = env_->template pool<SymEngine::Expression>().len(slot_);
            vec.resize(len);
            SymEngine::Expression* p =
                env_->template pool<SymEngine::Expression>().slot(slot_);
            for (unsigned i = 0; i < len; ++i) {
                if (!da::is_zero(p[i])) {
                    vec[i].init(vars, {p[i].get_basic()});
                } else {
                    vec[i].init(vars, {SymEngine::Expression(0).get_basic()});
                }
            }
        }
    }

    // eval_funs: produce a single LambdaRealDoubleVisitor for all elements
    void eval_funs(std::vector<SymEngine::RCP<const SymEngine::Basic>> vars,
                   SymEngine::LambdaRealDoubleVisitor& v) {
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            unsigned len = env_->template pool<SymEngine::Expression>().len(slot_);
            SymEngine::Expression* p =
                env_->template pool<SymEngine::Expression>().slot(slot_);
            std::vector<SymEngine::RCP<const SymEngine::Basic>> exprs;
            exprs.reserve(len);
            for (unsigned i = 0; i < len; ++i) {
                if (!da::is_zero(p[i])) {
                    exprs.push_back(p[i].get_basic());
                } else {
                    exprs.push_back(SymEngine::Expression(0).get_basic());
                }
            }
            v.init(vars, exprs);
        }
    }

#endif // DA_WITH_SYMBOLIC
};

// Static member definition
template<class T>
double DAVector<T>::eps = 1e-16;

// ===========================================================================
// Specialised element() for double — direct cast avoidance
// ===========================================================================

// ===========================================================================
// Free operators
// ===========================================================================

// --- DAVector op DAVector ---

template<class T>
DAVector<T> operator+(const DAVector<T>& a, const DAVector<T>& b) {
    check_env(a.env_, b.env_);
    DAVector<T> res;
    // res is alloc'd from current env; we need it in a's env
    // For single-env use (typical), current env == a.env_.
    // Force res to be in a's env:
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_copy(a.env_->layout(), a.env_->template pool<T>(), a.slot_, res.slot_);
    detail::ad_add(a.env_->layout(), a.env_->template pool<T>(), res.slot_, b.slot_);
    return res;
}

template<class T>
DAVector<T> operator-(const DAVector<T>& a, const DAVector<T>& b) {
    check_env(a.env_, b.env_);
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_copy(a.env_->layout(), a.env_->template pool<T>(), a.slot_, res.slot_);
    detail::ad_sub(a.env_->layout(), a.env_->template pool<T>(), res.slot_, b.slot_);
    return res;
}

template<class T>
DAVector<T> operator*(const DAVector<T>& a, const DAVector<T>& b) {
    check_env(a.env_, b.env_);
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_mult(a.env_->layout(), a.env_->template pool<T>(), a.slot_, b.slot_, res.slot_);
    return res;
}

template<class T>
DAVector<T> operator/(const DAVector<T>& a, const DAVector<T>& b) {
    check_env(a.env_, b.env_);
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_div(a.env_->layout(), a.env_->template pool<T>(), a.slot_, b.slot_, res.slot_);
    return res;
}

// --- DAVector op double ---

template<class T>
DAVector<T> operator+(const DAVector<T>& a, double x) {
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_copy(a.env_->layout(), a.env_->template pool<T>(), a.slot_, res.slot_);
    detail::ad_add_const(a.env_->layout(), a.env_->template pool<T>(), res.slot_,
                         static_cast<T>(x));
    return res;
}

template<class T>
DAVector<T> operator+(double x, const DAVector<T>& a) {
    return a + x;
}

template<class T>
DAVector<T> operator-(const DAVector<T>& a, double x) {
    return a + (-x);
}

template<class T>
DAVector<T> operator-(double x, const DAVector<T>& a) {
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    // res = -a + x
    detail::ad_mult_c(a.env_->layout(), a.env_->template pool<T>(), a.slot_,
                      static_cast<T>(-1.0), res.slot_);
    detail::ad_add_const(a.env_->layout(), a.env_->template pool<T>(), res.slot_,
                         static_cast<T>(x));
    return res;
}

template<class T>
DAVector<T> operator*(const DAVector<T>& a, double x) {
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_mult_c(a.env_->layout(), a.env_->template pool<T>(), a.slot_,
                      static_cast<T>(x), res.slot_);
    return res;
}

template<class T>
DAVector<T> operator*(double x, const DAVector<T>& a) {
    return a * x;
}

template<class T>
DAVector<T> operator/(const DAVector<T>& a, double x) {
    // Below DBL_MIN, x is zero or subnormal and 1/x overflows to inf.
    if (std::abs(x) < std::numeric_limits<double>::min())
        throw std::invalid_argument("da::operator/: divide by zero or a subnormal number");
    return a * (1.0 / x);
}

template<class T>
DAVector<T> operator/(double x, const DAVector<T>& a) {
    DAVector<T> res;
    if (res.env_ != a.env_) {
        res.env_->template pool<T>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<T>().alloc();
    }
    detail::ad_c_div(a.env_->layout(), a.env_->template pool<T>(), a.slot_,
                     static_cast<T>(x), res.slot_);
    return res;
}

// --- Unary ---

template<class T>
DAVector<T> operator+(const DAVector<T>& a) {
    return a;
}

template<class T>
DAVector<T> operator-(const DAVector<T>& a) {
    return a * (-1.0);
}

// --- Equality ---

template<class T>
bool operator==(const DAVector<T>& a, const DAVector<T>& b) {
    DAVector<T> c = a - b;
    return c.iszero();
}

// --- Stream output ---

template<class T>
std::ostream& operator<<(std::ostream& os, const DAVector<T>& v);

template<class T>
std::ostream& operator<<(std::ostream& os, const std::complex<DAVector<T>>& cd);

// ===========================================================================
// same_env free function for DAVector
// ===========================================================================
template<class T>
bool same_env(const DAVector<T>& a, const DAVector<T>& b) noexcept {
    return a.env_ == b.env_;
}


#ifdef DA_WITH_SYMBOLIC

// ==========================================================================
// SDA alias
// ==========================================================================
using SDA = DAVector<SymEngine::Expression>;

// ==========================================================================
// Free operators: DAVector<Expression> op Expression scalar (and reverse)
// ==========================================================================

inline SDA operator+(const SDA& a, SymEngine::Expression x) {
    SDA res;
    if (res.env_ != a.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_copy(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                    a.slot_, res.slot_);
    detail::ad_add_const(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                         res.slot_, x);
    return res;
}
inline SDA operator+(SymEngine::Expression x, const SDA& a) { return a + x; }

inline SDA operator-(const SDA& a, SymEngine::Expression x) {
    return a + (SymEngine::Expression(-1) * x);
}
inline SDA operator-(SymEngine::Expression x, const SDA& a) {
    SDA res;
    if (res.env_ != a.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_mult_c(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                      a.slot_, SymEngine::Expression(-1), res.slot_);
    detail::ad_add_const(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                         res.slot_, x);
    return res;
}

inline SDA operator*(const SDA& a, SymEngine::Expression x) {
    SDA res;
    if (res.env_ != a.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_mult_c(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                      a.slot_, x, res.slot_);
    return res;
}
inline SDA operator*(SymEngine::Expression x, const SDA& a) { return a * x; }

inline SDA operator/(const SDA& a, SymEngine::Expression x) {
    SDA res;
    if (res.env_ != a.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_copy(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                    a.slot_, res.slot_);
    detail::ad_div_c(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                     res.slot_, x);
    return res;
}
inline SDA operator/(SymEngine::Expression x, const SDA& a) {
    SDA res;
    if (res.env_ != a.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = a.env_;
        res.slot_ = a.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_c_div(a.env_->layout(), a.env_->template pool<SymEngine::Expression>(),
                     a.slot_, x, res.slot_);
    return res;
}

#endif // DA_WITH_SYMBOLIC

} // namespace da
