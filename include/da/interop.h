/**
 * @file interop.h
 * @brief NDA <-> SDA interoperability: promote, evaluate, mixed operators.
 *
 * @details Stage 7 per DEVELOPMENT_PLAN.md A.7.
 *   - promote(NDA) -> SDA: slot copy with expr[i] = Expression(dbl[i])
 *   - Mixed +,-,*,/ for (NDA,SDA) and (SDA,NDA) -> SDA
 *   - evaluate(SDA, map<symbol,double>) -> NDA
 *   - convert/import: double<->Expression, same-layout
 *
 * All declarations and definitions are guarded by DA_WITH_SYMBOLIC.
 * Include this header via the umbrella da/da.h.
 */

#pragma once

#ifdef DA_WITH_SYMBOLIC

#include "da/davector.h"
#include "da/engine.h"
#include "da/env.h"
#include "da/symbolic_ops.h"
#include <symengine/expression.h>
#include <symengine/eval_double.h>
#include <map>
#include <stdexcept>

namespace da {

// ===========================================================================
// promote: NDA -> SDA (same env, same layout, coefficient-by-coefficient)
// ===========================================================================

/**
 * @brief Promote a numerical DA vector to a symbolic DA vector.
 *
 * Allocates a new SDA slot in the SAME env as the source NDA,
 * copies each double coefficient as Expression(double_val),
 * and copies the slot length.
 *
 * The source NDA is unchanged.  The result SDA is in the same env.
 */
inline SDA promote(const NDA& src) {
    using E = SymEngine::Expression;
    DAEnv* env = src.env_;
    Pool<E>& epool = env->template pool<E>();
    Pool<double>& dpool = env->template pool<double>();

    SDA res;
    // Ensure res is in src's env (for multi-env safety)
    if (res.env_ != env) {
        res.env_->template pool<E>().free(res.slot_);
        res.env_ = env;
        res.slot_ = epool.alloc();
    }

    unsigned len = dpool.len(src.slot_);
    const double* src_ptr = dpool.slot(src.slot_);
    E* dst_ptr = epool.slot(res.slot_);

    // Zero-fill first (non-POD path: assign Expression(0))
    unsigned full_len = epool.full_len();
    for (unsigned i = 0; i < full_len; ++i)
        dst_ptr[i] = E(0);

    // Copy with promotion: use exact integer when possible to preserve
    // rational arithmetic (e.g., base vector coefficient 1.0 -> integer 1).
    for (unsigned i = 0; i < len; ++i) {
        double v = src_ptr[i];
        long iv = static_cast<long>(v);
        if (v == static_cast<double>(iv))
            dst_ptr[i] = E(iv);   // exact integer
        else
            dst_ptr[i] = E(v);    // float fallback
    }

    epool.set_len(res.slot_, len);
    return res;
}

// ===========================================================================
// Mixed operators: (NDA, SDA) and (SDA, NDA) -> SDA
// ===========================================================================

inline SDA operator+(const NDA& a, const SDA& b) {
    check_env(a.env_, b.env_);
    return promote(a) + b;
}
inline SDA operator+(const SDA& a, const NDA& b) {
    check_env(a.env_, b.env_);
    return a + promote(b);
}

inline SDA operator-(const NDA& a, const SDA& b) {
    check_env(a.env_, b.env_);
    return promote(a) - b;
}
inline SDA operator-(const SDA& a, const NDA& b) {
    check_env(a.env_, b.env_);
    return a - promote(b);
}

inline SDA operator*(const NDA& a, const SDA& b) {
    check_env(a.env_, b.env_);
    return promote(a) * b;
}
inline SDA operator*(const SDA& a, const NDA& b) {
    check_env(a.env_, b.env_);
    return a * promote(b);
}

inline SDA operator/(const NDA& a, const SDA& b) {
    check_env(a.env_, b.env_);
    return promote(a) / b;
}
inline SDA operator/(const SDA& a, const NDA& b) {
    check_env(a.env_, b.env_);
    return a / promote(b);
}

// ===========================================================================
// evaluate: SDA + symbol map -> NDA
// ===========================================================================

/**
 * @brief Evaluate an SDA at given symbol values, producing an NDA.
 *
 * Each symbolic coefficient is evaluated numerically using SymEngine's
 * eval_double.  The result NDA is allocated in the SAME env as the SDA.
 *
 * @param src     Symbolic DA vector to evaluate.
 * @param values  Map from SymEngine::RCP<Basic> (symbol) -> double.
 * @return NDA with the same length, each coefficient = numeric value.
 */
inline NDA evaluate(const SDA& src,
                    const SymEngine::map_basic_basic& values) {
    using E = SymEngine::Expression;
    DAEnv* env = src.env_;
    Pool<E>& epool = env->template pool<E>();
    Pool<double>& dpool = env->template pool<double>();

    NDA res;
    if (res.env_ != env) {
        res.env_->template pool<double>().free(res.slot_);
        res.env_ = env;
        res.slot_ = dpool.alloc();
    }

    unsigned len = epool.len(src.slot_);
    const E* src_ptr = epool.slot(src.slot_);
    double* dst_ptr = dpool.slot(res.slot_);

    // Zero fill
    unsigned full_len = dpool.full_len();
    for (unsigned i = 0; i < full_len; ++i)
        dst_ptr[i] = 0.0;

    // Evaluate each coefficient
    for (unsigned i = 0; i < len; ++i) {
        E ev = src_ptr[i].subs(values);
        dst_ptr[i] = SymEngine::eval_double(*ev.get_basic());
    }
    dpool.set_len(res.slot_, len);
    return res;
}

/**
 * @brief Evaluate using eval_funs (LambdaRealDoubleVisitor) for efficiency.
 *
 * @param src     Symbolic DA vector.
 * @param vars    Ordered list of symbolic variables (RCP<Basic>).
 * @param inputs  Numeric values for the variables, in the same order.
 * @return NDA with numerically evaluated coefficients.
 */
inline NDA evaluate(const SDA& src,
                    const std::vector<SymEngine::RCP<const SymEngine::Basic>>& vars,
                    const std::vector<double>& inputs) {
    using E = SymEngine::Expression;
    DAEnv* env = src.env_;
    Pool<E>& epool = env->template pool<E>();
    Pool<double>& dpool = env->template pool<double>();

    NDA res;
    if (res.env_ != env) {
        res.env_->template pool<double>().free(res.slot_);
        res.env_ = env;
        res.slot_ = dpool.alloc();
    }

    unsigned len = epool.len(src.slot_);
    const E* src_ptr = epool.slot(src.slot_);
    double* dst_ptr = dpool.slot(res.slot_);

    unsigned full_len = dpool.full_len();
    for (unsigned i = 0; i < full_len; ++i)
        dst_ptr[i] = 0.0;

    // Build one visitor per coefficient and evaluate
    for (unsigned i = 0; i < len; ++i) {
        if (da::is_zero(src_ptr[i])) {
            dst_ptr[i] = 0.0;
        } else {
            SymEngine::LambdaRealDoubleVisitor visitor;
            visitor.init(vars, {src_ptr[i].get_basic()});
            dst_ptr[i] = visitor.call(inputs);
        }
    }
    dpool.set_len(res.slot_, len);
    return res;
}

// ===========================================================================
// Mixed operators: Expression scalar op NDA -> SDA
// These allow expressions like: sx + sda[0]*sda[0] where sx is Expression
// and sda[0] is NDA (from da::base).
// ===========================================================================

inline SDA operator+(SymEngine::Expression x, const NDA& b) {
    SDA sb = promote(b);
    SDA res;
    if (res.env_ != sb.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sb.env_;
        res.slot_ = sb.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_copy(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                    sb.slot_, res.slot_);
    detail::ad_add_const(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                         res.slot_, x);
    return res;
}
inline SDA operator+(const NDA& a, SymEngine::Expression x) {
    return x + a;
}

inline SDA operator-(SymEngine::Expression x, const NDA& b) {
    SDA sb = promote(b);
    SDA res;
    if (res.env_ != sb.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sb.env_;
        res.slot_ = sb.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_mult_c(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                      sb.slot_, SymEngine::Expression(-1), res.slot_);
    detail::ad_add_const(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                         res.slot_, x);
    return res;
}
inline SDA operator-(const NDA& a, SymEngine::Expression x) {
    SDA sa = promote(a);
    SDA res;
    if (res.env_ != sa.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sa.env_;
        res.slot_ = sa.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_copy(sa.env_->layout(), sa.env_->template pool<SymEngine::Expression>(),
                    sa.slot_, res.slot_);
    detail::ad_add_const(sa.env_->layout(), sa.env_->template pool<SymEngine::Expression>(),
                         res.slot_, SymEngine::Expression(-1) * x);
    return res;
}

inline SDA operator*(SymEngine::Expression x, const NDA& b) {
    SDA sb = promote(b);
    SDA res;
    if (res.env_ != sb.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sb.env_;
        res.slot_ = sb.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_mult_c(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                      sb.slot_, x, res.slot_);
    return res;
}
inline SDA operator*(const NDA& a, SymEngine::Expression x) {
    return x * a;
}

inline SDA operator/(const NDA& a, SymEngine::Expression x) {
    SDA sa = promote(a);
    SDA res;
    if (res.env_ != sa.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sa.env_;
        res.slot_ = sa.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_copy(sa.env_->layout(), sa.env_->template pool<SymEngine::Expression>(),
                    sa.slot_, res.slot_);
    detail::ad_div_c(sa.env_->layout(), sa.env_->template pool<SymEngine::Expression>(),
                     res.slot_, x);
    return res;
}
inline SDA operator/(SymEngine::Expression x, const NDA& b) {
    SDA sb = promote(b);
    SDA res;
    if (res.env_ != sb.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = sb.env_;
        res.slot_ = sb.env_->template pool<SymEngine::Expression>().alloc();
    }
    detail::ad_c_div(sb.env_->layout(), sb.env_->template pool<SymEngine::Expression>(),
                     sb.slot_, x, res.slot_);
    return res;
}

// ===========================================================================
// promote: complex<NDA> -> complex<SDA> (promote both real and imaginary parts)
// ===========================================================================

/**
 * @brief Promote a numerical complex DA vector to a symbolic complex DA vector.
 *
 * Promotes both the real and imaginary parts via promote(NDA) -> SDA.
 */
inline std::complex<SDA> promote(const std::complex<NDA>& src) {
    const NDA& re = reinterpret_cast<const NDA(&)[2]>(src)[0];
    const NDA& im = reinterpret_cast<const NDA(&)[2]>(src)[1];
    return std::complex<SDA>(promote(re), promote(im));
}

// ===========================================================================
// evaluate: complex<SDA> + symbol map -> complex<NDA>
// ===========================================================================

/**
 * @brief Evaluate a symbolic complex DA vector at given symbol values, producing a complex<NDA>.
 *
 * Evaluates both real and imaginary parts via evaluate(SDA, map) -> NDA.
 *
 * @note This is why complex symbolic work belongs in CSDA rather than in
 *       SymEngine::I inside coefficients: each part must evaluate to a real
 *       double. An SDA carrying I makes SymEngine raise "Not Implemented"
 *       here. See the CSDA alias in da.h.
 */
inline std::complex<NDA> evaluate(const std::complex<SDA>& src,
                                   const SymEngine::map_basic_basic& values) {
    const SDA& re = reinterpret_cast<const SDA(&)[2]>(src)[0];
    const SDA& im = reinterpret_cast<const SDA(&)[2]>(src)[1];
    return std::complex<NDA>(evaluate(re, values), evaluate(im, values));
}

/**
 * @brief Evaluate a symbolic complex DA vector using an ordered variable/value list.
 */
inline std::complex<NDA> evaluate(const std::complex<SDA>& src,
                                   const std::vector<SymEngine::RCP<const SymEngine::Basic>>& vars,
                                   const std::vector<double>& inputs) {
    const SDA& re = reinterpret_cast<const SDA(&)[2]>(src)[0];
    const SDA& im = reinterpret_cast<const SDA(&)[2]>(src)[1];
    return std::complex<NDA>(evaluate(re, vars, inputs), evaluate(im, vars, inputs));
}

} // namespace da

#endif // DA_WITH_SYMBOLIC
