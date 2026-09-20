/**
 * @file functions.cpp
 * @brief Math functions for DAVector<double>.
 *
 * @details Ported from ref/tpsa/src/da.cc.  Algorithms are unchanged;
 *   only the API surface changes (namespace da, DAVector<double>/NDA,
 *   kernels take Layout+Pool instead of globals).
 */

#include "da/da.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace da {

// ===========================================================================
// sqrt
// ===========================================================================
NDA sqrt(const NDA& v) {
    NDA res(0.0);
    if (res.env_ != v.env_) {
        res.env_->pool<double>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->pool<double>().alloc();
    }
    double c = static_cast<double>(v.con());
    if (std::abs(c) < std::numeric_limits<double>::min()) {
        if (!const_cast<NDA&>(v).iszero()) {
            std::cout << "Warning: sqrt not defined because the constant part of the DA vector is zero." << std::endl;
            detail::ad_const(v.env_->layout(), v.env_->pool<double>(), res.slot_, (double)NAN);
        }
    } else {
        detail::ad_sqrt(v.env_->layout(), v.env_->pool<double>(), v.slot_, res.slot_);
    }
    return res;
}

// ===========================================================================
// exp
// ===========================================================================
NDA exp(const NDA& v) {
    NDA res;
    if (res.env_ != v.env_) {
        res.env_->pool<double>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->pool<double>().alloc();
    }
    detail::ad_exp(v.env_->layout(), v.env_->pool<double>(), v.slot_, res.slot_);
    return res;
}

// ===========================================================================
// log
// ===========================================================================
NDA log(const NDA& v) {
    NDA res;
    if (res.env_ != v.env_) {
        res.env_->pool<double>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->pool<double>().alloc();
    }
    double c = static_cast<double>(v.con());
    if (std::abs(c) < std::numeric_limits<double>::min()) {
        detail::ad_const(v.env_->layout(), v.env_->pool<double>(), res.slot_,
                         -std::numeric_limits<double>::infinity());
    } else {
        detail::ad_log(v.env_->layout(), v.env_->pool<double>(), v.slot_, res.slot_);
    }
    return res;
}

// ===========================================================================
// sin
// ===========================================================================
NDA sin(const NDA& v) {
    double cons = static_cast<double>(v.con());
    NDA pol = v - cons;
    NDA high_order = pol;
    int ord = NDA::order();
    std::vector<double> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = std::sin(cons);
    coefs[1] = std::cos(cons);
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] = -coefs[static_cast<std::size_t>(i-2)] / (i*(i-1));
    NDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// cos
// ===========================================================================
NDA cos(const NDA& v) {
    double cons = static_cast<double>(v.con());
    NDA pol = v - cons;
    NDA high_order = pol;
    int ord = NDA::order();
    std::vector<double> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = std::cos(cons);
    coefs[1] = -std::sin(cons);
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] = -coefs[static_cast<std::size_t>(i-2)] / (i*(i-1));
    NDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// tan
// ===========================================================================
NDA tan(const NDA& v) {
    return sin(v) / cos(v);
}

// ===========================================================================
// atan (port of da.cc atan — uses 4.4.34 / 4.4.42 from Abramowitz & Stegun)
// ===========================================================================
NDA atan(const NDA& v) {
    double da_cons = static_cast<double>(v.con());
    NDA da_high = (v - da_cons) / (v * da_cons + 1.0);
    NDA result(std::atan(da_cons));
    NDA da_high2 = da_high * da_high;
    for (double i = 1.0, j = -1.0; i < v.order() + 1; i += 2.0) {
        j *= -1.0;
        result += da_high * (j / i);
        da_high *= da_high2;
    }
    return result;
}

// ===========================================================================
// asin
// ===========================================================================
NDA asin(const NDA& v) {
    double cons = static_cast<double>(v.con());
    if (cons < -1.0 || cons > 1.0)
        throw std::domain_error("The constant part of the input DA vector of asin() should be within [-1,1].");
    NDA x = v / (1.0 + sqrt(1.0 - v * v));
    return 2.0 * atan(x);
}

// ===========================================================================
// acos
// ===========================================================================
NDA acos(const NDA& v) {
    double cons = static_cast<double>(v.con());
    if (cons < -1.0 || cons > 1.0)
        throw std::domain_error("The constant part of the input DA vector of acos() should be within [-1,1].");
    const double half_pi = 1.57079632679489661923132169163975144209858469968755291048;
    return NDA(half_pi) - asin(v);
}

// ===========================================================================
// sinh
// ===========================================================================
NDA sinh(const NDA& v) {
    double cons = static_cast<double>(v.con());
    NDA pol = v - cons;
    NDA high_order = pol;
    int ord = NDA::order();
    std::vector<double> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = std::sinh(cons);
    coefs[1] = std::cosh(cons);
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] = coefs[static_cast<std::size_t>(i-2)] / (i*(i-1));
    NDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// cosh
// ===========================================================================
NDA cosh(const NDA& v) {
    double cons = static_cast<double>(v.con());
    NDA pol = v - cons;
    NDA high_order = pol;
    int ord = NDA::order();
    std::vector<double> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = std::cosh(cons);
    coefs[1] = std::sinh(cons);
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] = coefs[static_cast<std::size_t>(i-2)] / (i*(i-1));
    NDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// tanh
// ===========================================================================
NDA tanh(const NDA& v) {
    return sinh(v) / cosh(v);
}

// ===========================================================================
// asinh
// ===========================================================================
NDA asinh(const NDA& v) {
    return log(v + sqrt(v * v + 1.0));
}

// ===========================================================================
// acosh
// ===========================================================================
NDA acosh(const NDA& v) {
    double cons = static_cast<double>(v.con());
    if (cons < 1.0)
        throw std::domain_error("Error in ACOSH: the constant part of the DA vector should be in [1, +inf).");
    return log(v + sqrt(v * v - 1.0));
}

// ===========================================================================
// atanh
// ===========================================================================
NDA atanh(const NDA& v) {
    double cons = static_cast<double>(v.con());
    if (cons <= -1.0 || cons >= 1.0)
        throw std::domain_error("Error in ATANH: the constant part of the DA vector should be in (-1,1).");
    return log((1.0 + v) / (1.0 - v)) / 2.0;
}

// ===========================================================================
// abs (DA version: max |coeff|)
// ===========================================================================
double abs(const NDA& v) {
    return v.norm();
}

double abs(const std::complex<NDA>& cv) {
    double r1 = abs(get_real(cv));
    double r2 = abs(get_imag(cv));
    return r1 > r2 ? r1 : r2;
}

// ===========================================================================
// pow helpers
// ===========================================================================

static NDA pow_pos(const NDA& v, int order) {
    NDA res;
    if (order == 0) {
        res.reset_const(1.0);
    } else if (order == 1) {
        res = v;
    } else if (order & 1) {
        res = v * pow_pos(v * v, order / 2);
    } else {
        res = pow_pos(v * v, order / 2);
    }
    return res;
}

NDA pow(const NDA& v, int order) {
    NDA res;
    if (order < 0) {
        double c = static_cast<double>(v.con());
        if (std::abs(c) < std::numeric_limits<double>::min()) {
            res.reset_const(std::numeric_limits<double>::infinity());
            return res;
        }
        res = pow_pos(v, -order);
        res = 1.0 / res;
    } else {
        res = pow_pos(v, order);
    }
    return res;
}

NDA pow(const NDA& v, double order) {
    if (std::floor(order) == order)
        return pow(v, static_cast<int>(order));
    double bas = static_cast<double>(v.con());
    NDA res;
    if (bas > 0) {
        res = exp(order * log(v));
    } else if (bas == 0.0) {
        res.reset_const(std::numeric_limits<double>::infinity());
    } else {
        res.reset_const(NAN);
    }
    return res;
}

// ===========================================================================
// erf
// ===========================================================================
NDA erf(const NDA& v) {
    const double coef = 1.1283791670955125585607; // 2/sqrt(pi)
    double cc = static_cast<double>(v.con());

    // Need to import da::base but we access via current env
    // Build base[0] manually: a DA var with constant cc
    NDA da0;
    {
        auto& env = da_current_env();
        da0.env_ = &env;
        detail::ad_var(env.layout(), env.pool<double>(), da0.slot_, 0.0, 0u);
    }

    NDA dal = NDA(cc) + da0;
    dal = exp(-1.0 * dal * dal);
    dal = da_int(dal, 0);
    NDA ada = v - NDA(cc);
    NDA da_erf_result;
    da_substitute(dal, 0, ada, da_erf_result);
    return da_erf_result * coef + NDA(std::erf(cc));
}

// ===========================================================================
// atan2
// ===========================================================================
NDA atan2(const NDA& y, const NDA& x) {
    NDA res;
    double cx = static_cast<double>(x.con());
    double cy = static_cast<double>(y.con());
    if (cx >= std::numeric_limits<double>::min()) {
        res = 2.0 * atan(y / (x + sqrt(x * x + y * y)));
    } else {
        if (std::abs(cy) > std::numeric_limits<double>::min()) {
            res = 2.0 * atan((sqrt(x * x + y * y) - x) / y);
        } else {
            if (std::abs(cx) < std::numeric_limits<double>::min()) {
                throw std::domain_error("Error: ATAN2 undefined.  Zero constant part for both the DA vectors.");
            } else {
                const double pi = 3.141592653589793238462643383279;
                res = atan(y / x) + NDA(pi);
            }
        }
    }
    return res;
}

// ===========================================================================
// atan2 — available for all DAVector<T>
// For T=double: use branch-based implementation (mirrors NDA atan2).
// For T=Expression: use the general-purpose formula that works symbolically.
// ===========================================================================
template<class T>
DAVector<T> atan2(const DAVector<T>& y, const DAVector<T>& x) {
    if constexpr (std::is_same<T, double>::value) {
        DAVector<T> res;
        double cx = static_cast<double>(x.con());
        double cy = static_cast<double>(y.con());
        if (cx >= std::numeric_limits<double>::min()) {
            res = 2.0 * atan(y / (x + sqrt(x * x + y * y)));
        } else {
            if (std::abs(cy) > std::numeric_limits<double>::min()) {
                res = 2.0 * atan((sqrt(x * x + y * y) - x) / y);
            } else {
                if (std::abs(cx) < std::numeric_limits<double>::min()) {
                    throw std::domain_error("Error: ATAN2 undefined. Zero constant part for both the DA vectors.");
                } else {
                    const double pi = 3.141592653589793238462643383279;
                    res = atan(y / x) + DAVector<T>(pi);
                }
            }
        }
        return res;
    } else {
        // Symbolic path: use the half-angle formula, valid when denominator is non-zero.
        // atan2(y, x) = 2 * atan(y / (sqrt(x^2 + y^2) + x))
        // This works as long as the constant part of (sqrt(x^2+y^2)+x) is non-zero.
        return 2.0 * atan(y / (x + sqrt(x * x + y * y)));
    }
}

// ===========================================================================
// Complex DA functions — template over DAVector<T>
// (ported from ref/tpsa/src/da.cc; body unchanged, NDA -> DAVector<T>)
// Explicitly instantiated below for T=double and (guarded) T=Expression.
// ===========================================================================

template<class T>
std::complex<DAVector<T>> exp(const std::complex<DAVector<T>>& c) {
    const DAVector<T>& rc = get_real(c);
    const DAVector<T>& ic = get_imag(c);
    const std::complex<double> ui(0.0, 1.0);
    return exp(rc) * (cos(ic) + ui * sin(ic));
}

template<class T>
std::complex<DAVector<T>> sqrt(const std::complex<DAVector<T>>& c) {
    const DAVector<T>& rc = get_real(c);
    const DAVector<T>& ic = get_imag(c);
    DAVector<T> r = sqrt(rc * rc + ic * ic);
    // Branch on whether constant part is on the negative real axis.
    // For T=double we can check numerically; for symbolic, use the general path.
    if constexpr (std::is_same<T, double>::value) {
        double ic_con = static_cast<double>(ic.con());
        double rc_con = static_cast<double>(rc.con());
        if (std::abs(ic_con) < std::numeric_limits<double>::min() &&
            rc_con < std::numeric_limits<double>::min()) {
            const std::complex<double> ui(0.0, 1.0);
            const double half_pi = 1.57079632679489661923132169163975144209858469968755291048;
            DAVector<T> theta(half_pi);
            if (!rc.iszero()) theta = theta + atan(ic / rc) / 2.0;
            return sqrt(r) * exp(ui * theta);
        }
    }
    // General (non-negative-real-axis) path — valid for both numeric and symbolic
    std::complex<DAVector<T>> cr = c + r;
    DAVector<T>& rcr = get_real(cr);
    DAVector<T>& icr = get_imag(cr);
    return sqrt(r) * cr / sqrt(rcr * rcr + icr * icr);
}

template<class T>
std::complex<DAVector<T>> log(const std::complex<DAVector<T>>& c) {
    std::complex<DAVector<T>> res;
    const DAVector<T>& rc = get_real(c);
    const DAVector<T>& ic = get_imag(c);
    // For T=double: handle the zero constant-part edge cases numerically.
    // For symbolic T: always use the general path (caller ensures non-zero).
    if constexpr (std::is_same<T, double>::value) {
        DAVector<T> r2 = rc * rc + ic * ic;
        double r2_con = static_cast<double>(r2.con());
        if (std::abs(r2_con) > std::numeric_limits<double>::min()) {
            get_real(res) = log(sqrt(rc * rc + ic * ic));
            get_imag(res) = atan2(ic, rc);
        } else {
            if (rc.iszero() && ic.iszero()) {
                get_real(res).reset_const(-std::numeric_limits<double>::infinity());
                get_imag(res) = DAVector<T>(0.0);
            } else {
                get_real(res) = DAVector<T>(0.0);
                get_imag(res) = atan2(ic, rc);
            }
        }
    } else {
        // Symbolic path: always compute log(|z|) + i*atan2(im, re)
        get_real(res) = log(sqrt(rc * rc + ic * ic));
        get_imag(res) = atan2(ic, rc);
    }
    return res;
}

template<class T>
std::complex<DAVector<T>> asin(const std::complex<DAVector<T>>& c) {
    const std::complex<double> ui(0.0, 1.0);
    return -ui * log(ui * c + sqrt(1.0 - c * c));
}

template<class T>
std::complex<DAVector<T>> acos(const std::complex<DAVector<T>>& c) {
    const std::complex<double> ui(0.0, 1.0);
    const double half_pi = 1.57079632679489661923132169163975144209858469968755291048;
    return DAVector<T>(half_pi) + ui * log(ui * c + sqrt(1.0 - c * c));
}

template<class T>
std::complex<DAVector<T>> atan(const std::complex<DAVector<T>>& c) {
    const std::complex<double> ui(0.0, 1.0);
    return ui * (log(1.0 - ui * c) - log(1.0 + ui * c)) / 2.0;
}

template<class T>
std::complex<DAVector<T>> asinh(const std::complex<DAVector<T>>& c) {
    return log(c + sqrt(c * c + 1.0));
}

template<class T>
std::complex<DAVector<T>> acosh(const std::complex<DAVector<T>>& c) {
    return log(c + sqrt(c * c - 1.0));
}

template<class T>
std::complex<DAVector<T>> atanh(const std::complex<DAVector<T>>& c) {
    return (log(1.0 + c) - log(1.0 - c)) / 2.0;
}

// Complex pow helpers — template
template<class T>
static std::complex<DAVector<T>> cpow_pos(const std::complex<DAVector<T>>& v, int order) {
    std::complex<DAVector<T>> res;
    if (order == 0) {
        get_real(res).reset_const(1.0);
        get_imag(res).reset_const(0.0);
    } else if (order == 1) {
        res = v;
    } else if (order & 1) {
        res = v * cpow_pos(v * v, order / 2);
    } else {
        res = cpow_pos(v * v, order / 2);
    }
    return res;
}

template<class T>
std::complex<DAVector<T>> pow(const std::complex<DAVector<T>>& v, int order) {
    if (order < 0) {
        std::complex<DAVector<T>> inv = 1.0 / v;
        return cpow_pos(inv, -order);
    }
    return cpow_pos(v, order);
}

template<class T>
std::complex<DAVector<T>> pow(const std::complex<DAVector<T>>& v, double order) {
    if (std::floor(order) == order)
        return pow(v, static_cast<int>(order));
    return exp(DAVector<T>(order) * log(v));
}

// ===========================================================================
// Explicit instantiations for T = double
// ===========================================================================
template std::complex<DAVector<double>> exp(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> sqrt(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> log(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> asin(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> acos(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> atan(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> asinh(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> acosh(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> atanh(const std::complex<DAVector<double>>&);
template std::complex<DAVector<double>> pow(const std::complex<DAVector<double>>&, int);
template std::complex<DAVector<double>> pow(const std::complex<DAVector<double>>&, double);


#ifdef DA_WITH_SYMBOLIC

// ===========================================================================
// Helper: allocate SDA result in v's env
// ===========================================================================
static SDA sda_result_in_env(const SDA& v) {
    SDA res;
    if (res.env_ != v.env_) {
        res.env_->template pool<SymEngine::Expression>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->template pool<SymEngine::Expression>().alloc();
    }
    return res;
}

// ===========================================================================
// sqrt (SDA)
// ===========================================================================
SDA sqrt(const SDA& v) {
    SDA res = sda_result_in_env(v);
    // Check constant term is non-zero symbolically
    const SymEngine::Expression& c0 = v.env_->template pool<SymEngine::Expression>().slot(v.slot_)[0];
    if (da::is_zero(c0)) {
        // warn and return NaN-constant SDA
        std::cerr << "Warning: sqrt not defined because constant part is zero." << std::endl;
        SymEngine::Expression nan_val(std::numeric_limits<double>::quiet_NaN());
        detail::ad_const(v.env_->layout(),
                         v.env_->template pool<SymEngine::Expression>(), res.slot_, nan_val);
    } else {
        detail::ad_sqrt(v.env_->layout(),
                        v.env_->template pool<SymEngine::Expression>(), v.slot_, res.slot_);
    }
    return res;
}

// ===========================================================================
// exp (SDA)
// ===========================================================================
SDA exp(const SDA& v) {
    SDA res = sda_result_in_env(v);
    detail::ad_exp(v.env_->layout(),
                   v.env_->template pool<SymEngine::Expression>(), v.slot_, res.slot_);
    return res;
}

// ===========================================================================
// log (SDA)
// ===========================================================================
SDA log(const SDA& v) {
    SDA res = sda_result_in_env(v);
    const SymEngine::Expression& c0 = v.env_->template pool<SymEngine::Expression>().slot(v.slot_)[0];
    if (da::is_zero(c0)) {
        SymEngine::Expression neg_inf(-std::numeric_limits<double>::infinity());
        detail::ad_const(v.env_->layout(),
                         v.env_->template pool<SymEngine::Expression>(), res.slot_, neg_inf);
    } else {
        detail::ad_log(v.env_->layout(),
                       v.env_->template pool<SymEngine::Expression>(), v.slot_, res.slot_);
    }
    return res;
}

// ===========================================================================
// sin (SDA) — mirrors NDA sin: Taylor expansion around constant part
// ===========================================================================
SDA sin(const SDA& v) {
    SymEngine::Expression cons = v.con();
    SDA pol = v - cons;
    SDA high_order = pol;
    int ord = SDA::order();
    using E = SymEngine::Expression;
    std::vector<E> coefs(static_cast<std::size_t>(ord + 1));
    // sin Taylor: coefs[0]=sin(c0), coefs[1]=cos(c0), coefs[k] = -coefs[k-2]/(k*(k-1))
    // Use exact integer divisor to preserve rational arithmetic in SymEngine.
    coefs[0] = E(SymEngine::sin(cons.get_basic()));
    coefs[1] = E(SymEngine::cos(cons.get_basic()));
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] =
            -coefs[static_cast<std::size_t>(i-2)] / E(i*(i-1));
    SDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// cos (SDA)
// ===========================================================================
SDA cos(const SDA& v) {
    SymEngine::Expression cons = v.con();
    SDA pol = v - cons;
    SDA high_order = pol;
    int ord = SDA::order();
    using E = SymEngine::Expression;
    std::vector<E> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = E(SymEngine::cos(cons.get_basic()));
    coefs[1] = E(SymEngine::Expression(-1) * SymEngine::sin(cons.get_basic()));
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] =
            -coefs[static_cast<std::size_t>(i-2)] / E(i*(i-1));
    SDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order *= pol;
    }
    return result;
}

// ===========================================================================
// tan (SDA)
// ===========================================================================
SDA tan(const SDA& v) {
    return sin(v) / cos(v);
}

// ===========================================================================
// atan (SDA)
// ===========================================================================
SDA atan(const SDA& v) {
    SymEngine::Expression da_cons = v.con();
    SDA da_high = (v - da_cons) / (v * da_cons + 1.0);
    SDA result(SymEngine::Expression(SymEngine::atan(da_cons.get_basic())));
    SDA da_high2 = da_high * da_high;
    for (double i = 1.0, j = -1.0; i < v.order() + 1; i += 2.0) {
        j *= -1.0;
        result += da_high * SymEngine::Expression(j / i);
        da_high *= da_high2;
    }
    return result;
}

// ===========================================================================
// asin (SDA)
// ===========================================================================
SDA asin(const SDA& v) {
    SDA x = v / (1.0 + sqrt(1.0 - v * v));
    return 2.0 * atan(x);
}

// ===========================================================================
// acos (SDA)
// ===========================================================================
SDA acos(const SDA& v) {
    const double half_pi = 1.57079632679489661923132169163975144209858469968755291048;
    return SDA(SymEngine::Expression(half_pi)) - asin(v);
}

// ===========================================================================
// sinh (SDA) — mirrors NDA direct Taylor expansion algorithm
// sinh(cons + pol) = sum_i coefs[i] * pol^i
// where coefs[0]=sinh(cons), coefs[1]=cosh(cons), coefs[i]=coefs[i-2]/(i*(i-1))
// Uses E(int) constructor to get exact SymEngine integer (preserves rational arithmetic).
// ===========================================================================
SDA sinh(const SDA& v) {
    using E = SymEngine::Expression;
    E cons = v.con();
    SDA pol = v - cons;          // purely non-constant part
    SDA high_order = pol;        // pol^1
    int ord = SDA::order();
    // symbolic Taylor coefficients — use exact integer divisors
    std::vector<E> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = E(SymEngine::sinh(cons.get_basic()));
    coefs[1] = E(SymEngine::cosh(cons.get_basic()));
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] =
            coefs[static_cast<std::size_t>(i-2)] / E(i*(i-1));
    SDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order = high_order * pol;
    }
    return result;
}

// ===========================================================================
// cosh (SDA) — mirrors NDA direct Taylor expansion algorithm
// cosh(cons + pol) = sum_i coefs[i] * pol^i
// where coefs[0]=cosh(cons), coefs[1]=sinh(cons), coefs[i]=coefs[i-2]/(i*(i-1))
// Uses E(int) constructor to get exact SymEngine integer (preserves rational arithmetic).
// ===========================================================================
SDA cosh(const SDA& v) {
    using E = SymEngine::Expression;
    E cons = v.con();
    SDA pol = v - cons;
    SDA high_order = pol;
    int ord = SDA::order();
    std::vector<E> coefs(static_cast<std::size_t>(ord + 1));
    coefs[0] = E(SymEngine::cosh(cons.get_basic()));
    coefs[1] = E(SymEngine::sinh(cons.get_basic()));
    for (int i = 2; i < ord + 1; ++i)
        coefs[static_cast<std::size_t>(i)] =
            coefs[static_cast<std::size_t>(i-2)] / E(i*(i-1));
    SDA result(coefs[0]);
    for (int i = 1; i < ord + 1; ++i) {
        result += coefs[static_cast<std::size_t>(i)] * high_order;
        high_order = high_order * pol;
    }
    return result;
}

// ===========================================================================
// tanh (SDA)
// ===========================================================================
SDA tanh(const SDA& v) {
    return sinh(v) / cosh(v);
}

// ===========================================================================
// pow (SDA, int)
// ===========================================================================
static SDA sda_pow_pos(const SDA& v, int order) {
    if (order == 0) {
        SDA res = sda_result_in_env(v);
        detail::ad_const(v.env_->layout(),
                         v.env_->template pool<SymEngine::Expression>(),
                         res.slot_, SymEngine::Expression(1));
        return res;
    }
    if (order == 1) return v;
    if (order & 1) return v * sda_pow_pos(v * v, order / 2);
    return sda_pow_pos(v * v, order / 2);
}

SDA pow(const SDA& v, int order) {
    if (order < 0) {
        SDA inv = SymEngine::Expression(1) / v;
        return sda_pow_pos(inv, -order);
    }
    return sda_pow_pos(v, order);
}

// ===========================================================================
// pow (SDA, double)
// ===========================================================================
SDA pow(const SDA& v, double order) {
    if (std::floor(order) == order) return pow(v, static_cast<int>(order));
    return exp(SDA(SymEngine::Expression(order)) * log(v));
}

// ===========================================================================
// erf (SDA) — mirrors the reference sda.cc implementation
// Uses: erf(cc + ada) = coef * integral(exp(-t^2)) substituting ada for t,
// plus erf(cc) as constant term.
// ===========================================================================
SDA erf(const SDA& v) {
    const double coef = 1.1283791670955125585607; // 2/sqrt(pi)
    using E = SymEngine::Expression;
    E cons = v.con();

    // Build SDA base variable for variable 0: constant cons, linear term 1
    // This mirrors da[0] in the reference (the 0-th base vector)
    // We build: dal = cons + base[0]  where base[0] is the first NDA base vector
    SDA da0 = cons + da::base[0];   // Expression + NDA -> SDA via interop.h
    SDA dal = exp(-1.0 * da0 * da0);
    dal = da_int(dal, 0u);
    SDA ada = v - cons;             // SDA - Expression
    SDA da_erf_result;
    da_substitute(dal, 0u, ada, da_erf_result);
    return da_erf_result * coef + SDA(E(SymEngine::erf(cons.get_basic())));
}

// ===========================================================================
// abs(complex<SDA>) — symbolic magnitude: sqrt(re^2 + im^2) as SDA
// ===========================================================================
SDA abs(const std::complex<SDA>& cv) {
    const SDA& re = get_real(cv);
    const SDA& im = get_imag(cv);
    return sqrt(re * re + im * im);
}

// ===========================================================================
// Explicit instantiations for T = Expression (symbolic complex DA functions)
// ===========================================================================
template std::complex<DAVector<SymEngine::Expression>> exp(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> sqrt(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> log(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> asin(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> acos(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> atan(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> asinh(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> acosh(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> atanh(const std::complex<DAVector<SymEngine::Expression>>&);
template std::complex<DAVector<SymEngine::Expression>> pow(const std::complex<DAVector<SymEngine::Expression>>&, int);
template std::complex<DAVector<SymEngine::Expression>> pow(const std::complex<DAVector<SymEngine::Expression>>&, double);

#endif // DA_WITH_SYMBOLIC

} // namespace da
