/**
 * @file da.h
 * @brief Umbrella header for the da library (numerical complete after Stage 5).
 *
 * @details Include this single header for all DA functionality.
 *   Provides: Layout, Pool<T>, DAEnv, DAVector<T>/NDA, math functions,
 *   da::base[] bases, and top-level da_* free functions.
 */

#pragma once

#include <cassert>

#include "da/monomial_scheme.h"
#include "da/layout.h"
#include "da/pool.h"
#include "da/env.h"
#include "da/engine.h"
#include "da/davector.h"
#include "da/functions.h"
#include "da/base.h"
#ifdef DA_WITH_SYMBOLIC
#  include "da/interop.h"
#endif

#include <complex>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

namespace da {

// ===========================================================================
// Type trait to constrain complex DAVector templates to known coefficient types.
// This prevents ambiguity with SymEngine::Expression's own operator overloads.
// ===========================================================================

template<class T> struct is_da_coeff : std::false_type {};
template<> struct is_da_coeff<double> : std::true_type {};
#ifdef DA_WITH_SYMBOLIC
template<> struct is_da_coeff<SymEngine::Expression> : std::true_type {};
#endif

// ===========================================================================
// Complex DA helpers — template over DAVector<T>
// (mirrors ref/tpsa/include/da.h; body unchanged, only NDA -> DAVector<T>)
// ===========================================================================

template<class T>
inline DAVector<T>& get_real(std::complex<DAVector<T>>& v) {
    return reinterpret_cast<DAVector<T>(&)[2]>(v)[0];
}
template<class T>
inline DAVector<T>& get_imag(std::complex<DAVector<T>>& v) {
    return reinterpret_cast<DAVector<T>(&)[2]>(v)[1];
}
template<class T>
inline const DAVector<T>& get_real(const std::complex<DAVector<T>>& v) {
    return reinterpret_cast<const DAVector<T>(&)[2]>(v)[0];
}
template<class T>
inline const DAVector<T>& get_imag(const std::complex<DAVector<T>>& v) {
    return reinterpret_cast<const DAVector<T>(&)[2]>(v)[1];
}

// ===========================================================================
// cd_copy helpers — template over DAVector<T>
// ===========================================================================
template<class T>
inline void cd_copy(std::complex<DAVector<T>>& vs, std::complex<DAVector<T>>& vo) {
    get_real(vo) = get_real(vs);
    get_imag(vo) = get_imag(vs);
}
template<class T>
inline void cd_copy(std::complex<double> vs, std::complex<DAVector<T>>& vo) {
    get_real(vo) = DAVector<T>(vs.real());
    get_imag(vo) = DAVector<T>(vs.imag());
}
template<class T>
inline void cd_copy(double x, std::complex<DAVector<T>>& vo) {
    get_real(vo) = DAVector<T>(x);
    get_imag(vo) = DAVector<T>(0.0);
}
// cd_copy from a scalar of type T (e.g., Expression constant) — symbolic-only path
template<class T, typename std::enable_if<!std::is_same<T, double>::value, int>::type = 0>
inline void cd_copy(const T& x, std::complex<DAVector<T>>& vo) {
    get_real(vo) = DAVector<T>(x);
    get_imag(vo) = DAVector<T>(T(0.0));
}

// ===========================================================================
// Complex<DAVector<T>> operators — template over T
// (ported from ref/tpsa/src/da.cc; body unchanged, NDA -> DAVector<T>)
// ===========================================================================

template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(const DAVector<T>& v, std::complex<double> c) {
    return std::complex<DAVector<T>>(c.real() + v, DAVector<T>(c.imag()));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(std::complex<double> c, const DAVector<T>& v) {
    return std::complex<DAVector<T>>(c.real() + v, DAVector<T>(c.imag()));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(const DAVector<T>& v, std::complex<double> c) {
    return std::complex<DAVector<T>>(v - c.real(), DAVector<T>(c.imag()));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(std::complex<double> c, const DAVector<T>& v) {
    return std::complex<DAVector<T>>(c.real() - v, DAVector<T>(c.imag()));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(const DAVector<T>& v, std::complex<double> c) {
    return std::complex<DAVector<T>>(v * c.real(), v * c.imag());
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(std::complex<double> c, const DAVector<T>& v) {
    return std::complex<DAVector<T>>(c.real() * v, c.imag() * v);
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(const DAVector<T>& v, std::complex<double> c) {
    std::complex<double> inv = 1.0 / c;
    return std::complex<DAVector<T>>(v * inv.real(), v * inv.imag());
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(std::complex<double> c, const DAVector<T>& v) {
    DAVector<T> d = 1.0 / v;
    return std::complex<DAVector<T>>(c.real() * d, c.imag() * d);
}

// complex<DAVector<T>> +-*/ double
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(const std::complex<DAVector<T>>& cd, double n) {
    return std::complex<DAVector<T>>(get_real(cd) + n, get_imag(cd));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(double n, const std::complex<DAVector<T>>& cd) {
    return cd + n;
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(const std::complex<DAVector<T>>& cd, double n) {
    return std::complex<DAVector<T>>(get_real(cd) - n, get_imag(cd));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(double n, const std::complex<DAVector<T>>& cd) {
    return std::complex<DAVector<T>>(n - get_real(cd), -1.0 * get_imag(cd));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(const std::complex<DAVector<T>>& cd, double n) {
    return std::complex<DAVector<T>>(get_real(cd) * n, get_imag(cd) * n);
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(double n, const std::complex<DAVector<T>>& cd) {
    return cd * n;
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(const std::complex<DAVector<T>>& cd, double n) {
    return std::complex<DAVector<T>>(get_real(cd) / n, get_imag(cd) / n);
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(double n, const std::complex<DAVector<T>>& cd) {
    DAVector<T> rv = get_real(cd);
    DAVector<T> iv = get_imag(cd);
    DAVector<T> v = 1.0 / (rv * rv + iv * iv);
    return std::complex<DAVector<T>>(n * v * rv, -n * v * iv);
}

// complex<DAVector<T>> +-*/ complex<double>
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(const std::complex<DAVector<T>>& cd, std::complex<double> c) {
    return std::complex<DAVector<T>>(get_real(cd) + c.real(), get_imag(cd) + c.imag());
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(std::complex<double> c, const std::complex<DAVector<T>>& cd) {
    return cd + c;
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(const std::complex<DAVector<T>>& cd, std::complex<double> c) {
    return std::complex<DAVector<T>>(get_real(cd) - c.real(), get_imag(cd) - c.imag());
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(std::complex<double> c, const std::complex<DAVector<T>>& cd) {
    return std::complex<DAVector<T>>(c.real() - get_real(cd), c.imag() - get_imag(cd));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(const std::complex<DAVector<T>>& cd, std::complex<double> c) {
    auto vr = cd * c.real();
    auto vi = cd * c.imag();
    return std::complex<DAVector<T>>(get_real(vr) - get_imag(vi), get_imag(vr) + get_real(vi));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(std::complex<double> c, const std::complex<DAVector<T>>& cd) {
    return cd * c;
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(const std::complex<DAVector<T>>& cd, std::complex<double> c) {
    return (1.0 / c) * cd;
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(std::complex<double> c, const std::complex<DAVector<T>>& cd) {
    std::complex<DAVector<T>> inv = 1.0 / cd;
    return c * inv;
}

// complex<DAVector<T>> +-*/ complex<DAVector<T>>
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator+(const std::complex<DAVector<T>>& a, const std::complex<DAVector<T>>& b) {
    return std::complex<DAVector<T>>(get_real(a) + get_real(b), get_imag(a) + get_imag(b));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator-(const std::complex<DAVector<T>>& a, const std::complex<DAVector<T>>& b) {
    return std::complex<DAVector<T>>(get_real(a) - get_real(b), get_imag(a) - get_imag(b));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator*(const std::complex<DAVector<T>>& a, const std::complex<DAVector<T>>& b) {
    return std::complex<DAVector<T>>(
        get_real(a) * get_real(b) - get_imag(a) * get_imag(b),
        get_real(a) * get_imag(b) + get_imag(a) * get_real(b));
}
template<class T, typename std::enable_if<is_da_coeff<T>::value, int>::type = 0>
inline std::complex<DAVector<T>> operator/(const std::complex<DAVector<T>>& a, const std::complex<DAVector<T>>& b) {
    return a * (1.0 / b);
}

// ===========================================================================
// Type aliases for complex DA vectors
// ===========================================================================
using CNDA = std::complex<NDA>;
#ifdef DA_WITH_SYMBOLIC
/**
 * @brief Symbolic complex DA vector: real+imaginary pair of SDA.
 *
 * This is the supported way to do complex arithmetic with symbolic DA.
 * The imaginary unit belongs in this std::complex wrapper, never inside a
 * coefficient: an SDA holding SymEngine::I compiles and looks plausible, but
 * the real and imaginary parts can no longer be separated, evaluate() cannot
 * convert it back to an NDA (SymEngine raises "Not Implemented", since NDA
 * coefficients are double), and mixing the two representations double-counts
 * the imaginary unit.
 *
 * @note The complex function set matches CNDA and ref/tpsa exactly:
 *       exp, sqrt, log, asin, acos, atan, asinh, acosh, atanh, pow, abs.
 *       There are deliberately no complex sin/cos/tan/sinh/cosh/tanh
 *       overloads; the scalar SDA versions do exist.
 */
using CSDA = std::complex<SDA>;

// CSDA (op) Expression: a symbolic scalar is real, so it acts on the real part
// for + and -, and on both parts for * and /.
inline CSDA operator+(const CSDA& z, const SymEngine::Expression& e) { return CSDA(get_real(z) + e, get_imag(z)); }
inline CSDA operator+(const SymEngine::Expression& e, const CSDA& z) { return z + e; }
inline CSDA operator-(const CSDA& z, const SymEngine::Expression& e) { return CSDA(get_real(z) - e, get_imag(z)); }
inline CSDA operator-(const SymEngine::Expression& e, const CSDA& z) { return CSDA(e - get_real(z), -get_imag(z)); }
inline CSDA operator*(const CSDA& z, const SymEngine::Expression& e) { return CSDA(get_real(z) * e, get_imag(z) * e); }
inline CSDA operator*(const SymEngine::Expression& e, const CSDA& z) { return z * e; }
inline CSDA operator/(const CSDA& z, const SymEngine::Expression& e) { return CSDA(get_real(z) / e, get_imag(z) / e); }
inline CSDA operator/(const SymEngine::Expression& e, const CSDA& z) { return e * (1.0 / z); }
#endif

// ===========================================================================
// Top-level free functions — delegate to current env (mirrors ref/tpsa/include/da.h)
// ===========================================================================

/// Initialize the default DA environment.
// (already declared in env.h; re-declared with default for Base setup)

/// Set the cutoff value for DA coefficients.
inline void da_set_eps(double eps) {
    if (eps > 0) NDA::eps = eps;
    else std::cout << "Warning: negative value of eps is ignored in da_set_eps!" << std::endl;
}

/// Number of DA vectors currently allocated.
inline int da_count() {
    return static_cast<int>(da_current_env().pool<double>().count());
}

/// Number of free slots remaining in the pool.
inline int da_remain() {
    return static_cast<int>(da_current_env().pool<double>().remain());
}

/// Total pool size.
inline int da_poolsize() {
    return static_cast<int>(da_current_env().pool<double>().poolsize());
}

/// Full length of a DA vector.
inline int da_full_length() {
    return static_cast<int>(da_current_env().layout().full_len());
}

/// Return the exponent vector (orders) of the i-th element.
inline const std::vector<int>& da_element_orders(int i) {
    return da_current_env().layout().orders_ref(static_cast<unsigned int>(i));
}

/// Temporarily lower the DA order.
inline int da_change_order(unsigned int new_order) {
    int r = da_current_env().layout().change_order(new_order);
    if (r == 1)
        std::cout << "WARNING: Given order is too large. DA order does not change!" << std::endl;
    return r;
}

/// Restore the original DA order.
inline int da_restore_order() {
    da_current_env().layout().restore_order();
    return 0;
}

/// Take derivative of a DA vector w.r.t. base base_id.
inline void da_der(const NDA& v, unsigned int base_id, NDA& result) {
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base out of limits in da_der!");
    detail::ad_der(v.env_->layout(), v.env_->pool<double>(), v.slot_, base_id, result.slot_);
}

inline NDA da_der(const NDA& v, unsigned int base_id) {
    NDA res;
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base out of limits in da_der!");
    detail::ad_der(v.env_->layout(), v.env_->pool<double>(), v.slot_, base_id, res.slot_);
    return res;
}

/// Integrate a DA vector w.r.t. base base_id.
inline void da_int(const NDA& v, unsigned int base_id, NDA& result) {
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base out of limits in da_int!");
    detail::ad_int(v.env_->layout(), v.env_->pool<double>(), v.slot_, base_id, result.slot_);
}

inline NDA da_int(const NDA& v, unsigned int base_id) {
    NDA res;
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base out of limits in da_int!");
    detail::ad_int(v.env_->layout(), v.env_->pool<double>(), v.slot_, base_id, res.slot_);
    return res;
}

/// Substitute a constant for a base in a DA vector.
inline void da_substitute_const(const NDA& iv, unsigned int base_id, double x, NDA& ov) {
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base id out of range in da_substitute_const!");
    // Use single-vector substitute with a constant slot
    NDA tmp(x);
    detail::ad_substitute(iv.env_->layout(), iv.env_->pool<double>(),
                          iv.slot_, base_id, tmp.slot_, ov.slot_);
}

/// Substitute a DA vector for a base.
inline void da_substitute(const NDA& iv, unsigned int base_id, const NDA& v, NDA& ov) {
    assert(base_id < static_cast<unsigned>(NDA::dim()) && "Base id out of range in da_substitute!");
    detail::ad_substitute(iv.env_->layout(), iv.env_->pool<double>(),
                          iv.slot_, base_id, v.slot_, ov.slot_);
}

/// Substitute a constant for a base.
inline void da_substitute(const NDA& iv, unsigned int base_id, double x, NDA& ov) {
    da_substitute_const(iv, base_id, x, ov);
}

/// Substitute multiple DA vectors for multiple bases.
inline void da_substitute(const NDA& iv,
                          std::vector<unsigned int>& base_id,
                          std::vector<NDA>& v,
                          NDA& ov) {
    std::set<unsigned int> check(base_id.begin(), base_id.end());
    assert(check.size() == base_id.size() && "Duplicate base id in da_substitute!");
    for (auto chk_id : check)
        if (chk_id >= static_cast<unsigned>(NDA::dim()))
            throw std::invalid_argument("Base id out of range in da_substitute!");
    std::vector<unsigned> ad_v, ad_base_id;
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto id : base_id) ad_base_id.push_back(id);
    std::vector<unsigned> ad_iv = { iv.slot_ };
    std::vector<unsigned> ad_ov = { ov.slot_ };
    detail::ad_substitute(iv.env_->layout(), iv.env_->pool<double>(),
                          ad_iv, ad_base_id, ad_v, ad_ov);
}

/// Substitute multiple DA vectors for multiple bases in a group of vectors.
inline void da_substitute(std::vector<NDA>& ivecs,
                          std::vector<unsigned int>& base_id,
                          std::vector<NDA>& v,
                          std::vector<NDA>& ovecs) {
    std::set<unsigned int> check(base_id.begin(), base_id.end());
    assert(check.size() == base_id.size() && "Duplicate base id in da_substitute!");
    for (auto chk_id : check)
        if (chk_id >= static_cast<unsigned>(NDA::dim()))
            throw std::invalid_argument("Base id out of range in da_substitute!");
    std::vector<unsigned> ad_iv, ad_v, ad_ov, ad_base_id;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto& vv : ovecs) ad_ov.push_back(vv.slot_);
    for (auto id : base_id) ad_base_id.push_back(id);
    if (ivecs.empty()) return;
    detail::ad_substitute(ivecs[0].env_->layout(), ivecs[0].env_->pool<double>(),
                          ad_iv, ad_base_id, ad_v, ad_ov);
}

/// Composition of DA vectors with DA vectors.
inline void da_composition(std::vector<NDA>& ivecs, std::vector<NDA>& v, std::vector<NDA>& ovecs) {
    std::vector<unsigned> ad_iv, ad_v, ad_ov;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto& vv : ovecs) ad_ov.push_back(vv.slot_);
    if (ivecs.empty()) return;
    detail::ad_composition(ivecs[0].env_->layout(), ivecs[0].env_->pool<double>(),
                           ad_iv, ad_v, ad_ov);
}

/// Composition of DA vectors with double values.
inline void da_composition(std::vector<NDA>& ivecs, std::vector<double>& v, std::vector<double>& ovecs) {
    std::vector<unsigned> ad_iv;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    if (ivecs.empty()) return;
    detail::ad_composition(ivecs[0].env_->layout(), ivecs[0].env_->pool<double>(),
                           ad_iv, v, ovecs);
}

/// Composition of DA vectors with complex double values.
inline void da_composition(std::vector<NDA>& ivecs,
                           std::vector<std::complex<double>>& v,
                           std::vector<std::complex<double>>& ovecs) {
    std::vector<unsigned> ad_iv;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    if (ivecs.empty()) return;
    detail::ad_composition(ivecs[0].env_->layout(), ivecs[0].env_->pool<double>(),
                           ad_iv, v, ovecs);
}

// ===========================================================================
// cd_composition (complex DA composition — ported from da.cc as-is)
// Template over DAVector<T>; explicitly instantiated for T=double (always)
// and T=Expression (guarded). The NDA-specific overloads remain as before.
// ===========================================================================

// Forward declaration of complex_da_pow_int_pos (defined in base.cpp)
template<class T>
std::complex<DAVector<T>>& complex_da_pow_int_pos(const std::complex<DAVector<T>>& iv,
                                                   std::vector<std::complex<DAVector<T>>>& power_v,
                                                   int order, int order_rec,
                                                   std::vector<bool>* inited = nullptr);

template<class T>
void cd_composition(std::vector<DAVector<T>>& ivecs,
                    std::vector<std::complex<DAVector<T>>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs);

template<class T>
void cd_composition(std::vector<std::complex<DAVector<T>>>& ivecs,
                    std::vector<std::complex<DAVector<T>>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs);

template<class T>
void cd_composition(std::vector<std::complex<DAVector<T>>>& ivecs,
                    std::vector<DAVector<T>>& v,
                    std::vector<std::complex<DAVector<T>>>& ovecs);

// ===========================================================================
// inv_map
// ===========================================================================
void inv_map(std::vector<NDA>& ivecs, int dim, std::vector<NDA>& ovecs);

// ===========================================================================
// File I/O helpers (ported from da.cc)
// ===========================================================================
std::string trim_whitespace(std::string input_line);
bool read_da_from_file(std::string filename, NDA& d);
bool read_cd_from_file(std::string filename, std::complex<NDA>& cd);
NDA devide_by_element(NDA& t, NDA& b);
bool compare_da_vectors(NDA& a, NDA& b, double eps = 1e-15);
bool compare_da_with_file(std::string filename, NDA& d, double eps = 1e-15);
bool compare_cd_vectors(std::complex<NDA>& a, std::complex<NDA>& b, double eps = 1e-15);
bool compare_cd_with_file(std::string filename, std::complex<NDA>& d, double eps = 1e-15);

#ifdef DA_WITH_SYMBOLIC
// ===========================================================================
// SDA top-level free functions (symbolic DA)
// Mirrors the NDA versions above, adapted for SDA = DAVector<Expression>.
// ===========================================================================

/// Take derivative of an SDA vector w.r.t. base base_id.
inline void da_der(const SDA& v, unsigned int base_id, SDA& result) {
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base out of limits in da_der (SDA)!");
    detail::ad_der(v.env_->layout(), v.env_->template pool<E>(), v.slot_, base_id, result.slot_);
}

inline SDA da_der(const SDA& v, unsigned int base_id) {
    SDA res;
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base out of limits in da_der (SDA)!");
    if (res.env_ != v.env_) {
        res.env_->template pool<E>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->template pool<E>().alloc();
    }
    detail::ad_der(v.env_->layout(), v.env_->template pool<E>(), v.slot_, base_id, res.slot_);
    return res;
}

/// Integrate an SDA vector w.r.t. base base_id.
inline void da_int(const SDA& v, unsigned int base_id, SDA& result) {
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base out of limits in da_int (SDA)!");
    detail::ad_int(v.env_->layout(), v.env_->template pool<E>(), v.slot_, base_id, result.slot_);
}

inline SDA da_int(const SDA& v, unsigned int base_id) {
    SDA res;
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base out of limits in da_int (SDA)!");
    if (res.env_ != v.env_) {
        res.env_->template pool<E>().free(res.slot_);
        res.env_ = v.env_;
        res.slot_ = v.env_->template pool<E>().alloc();
    }
    detail::ad_int(v.env_->layout(), v.env_->template pool<E>(), v.slot_, base_id, res.slot_);
    return res;
}

/// Substitute a constant for a base in an SDA vector.
inline void da_substitute_const(const SDA& iv, unsigned int base_id, double x, SDA& ov) {
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base id out of range in da_substitute_const (SDA)!");
    SDA tmp{E(x)};
    detail::ad_substitute(iv.env_->layout(), iv.env_->template pool<E>(),
                          iv.slot_, base_id, tmp.slot_, ov.slot_);
}

/// Substitute an SDA vector for a base.
inline void da_substitute(const SDA& iv, unsigned int base_id, const SDA& v, SDA& ov) {
    using E = SymEngine::Expression;
    assert(base_id < static_cast<unsigned>(SDA::dim()) && "Base id out of range in da_substitute (SDA)!");
    detail::ad_substitute(iv.env_->layout(), iv.env_->template pool<E>(),
                          iv.slot_, base_id, v.slot_, ov.slot_);
}

/// Substitute multiple SDA vectors for multiple bases.
inline void da_substitute(const SDA& iv,
                          std::vector<unsigned int>& base_id,
                          std::vector<SDA>& v,
                          SDA& ov) {
    using E = SymEngine::Expression;
    std::set<unsigned int> check(base_id.begin(), base_id.end());
    assert(check.size() == base_id.size() && "Duplicate base id in da_substitute (SDA)!");
    for (auto chk_id : check)
        if (chk_id >= static_cast<unsigned>(SDA::dim()))
            throw std::invalid_argument("Base id out of range in da_substitute (SDA)!");
    std::vector<unsigned> ad_v, ad_base_id;
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto id : base_id) ad_base_id.push_back(id);
    std::vector<unsigned> ad_iv = { iv.slot_ };
    std::vector<unsigned> ad_ov = { ov.slot_ };
    detail::ad_substitute(iv.env_->layout(), iv.env_->template pool<E>(),
                          ad_iv, ad_base_id, ad_v, ad_ov);
}

/// Substitute multiple SDA vectors for multiple bases in a group of SDA vectors.
inline void da_substitute(std::vector<SDA>& ivecs,
                          std::vector<unsigned int>& base_id,
                          std::vector<SDA>& v,
                          std::vector<SDA>& ovecs) {
    using E = SymEngine::Expression;
    std::set<unsigned int> check(base_id.begin(), base_id.end());
    assert(check.size() == base_id.size() && "Duplicate base id in da_substitute (SDA)!");
    for (auto chk_id : check)
        if (chk_id >= static_cast<unsigned>(SDA::dim()))
            throw std::invalid_argument("Base id out of range in da_substitute (SDA)!");
    std::vector<unsigned> ad_iv, ad_v, ad_ov, ad_base_id;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto& vv : ovecs) ad_ov.push_back(vv.slot_);
    for (auto id : base_id) ad_base_id.push_back(id);
    if (ivecs.empty()) return;
    detail::ad_substitute(ivecs[0].env_->layout(), ivecs[0].env_->template pool<E>(),
                          ad_iv, ad_base_id, ad_v, ad_ov);
}

/// Composition of SDA vectors with SDA vectors.
inline void da_composition(std::vector<SDA>& ivecs, std::vector<SDA>& v, std::vector<SDA>& ovecs) {
    using E = SymEngine::Expression;
    std::vector<unsigned> ad_iv, ad_v, ad_ov;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    for (auto& vv : v) ad_v.push_back(vv.slot_);
    for (auto& vv : ovecs) ad_ov.push_back(vv.slot_);
    if (ivecs.empty()) return;
    detail::ad_composition(ivecs[0].env_->layout(), ivecs[0].env_->template pool<E>(),
                           ad_iv, ad_v, ad_ov);
}

/// Composition of SDA vectors with double values -> Expression scalar results.
/// (Mirrors ref/tpsa_sym/include/sda.h: da_composition(vector<SDA>, vector<double>, vector<Expression>))
inline void da_composition(std::vector<SDA>& ivecs,
                           std::vector<double>& v,
                           std::vector<SymEngine::Expression>& ovecs) {
    using E = SymEngine::Expression;
    if (ivecs.empty()) return;
    DAEnv* env = ivecs[0].env_;
    unsigned nvars = static_cast<unsigned>(v.size());
    // Build temporary SDA constants for each double input value
    std::vector<SDA> tmp_v(nvars), tmp_ov(ivecs.size());
    for (unsigned i = 0; i < nvars; ++i) {
        tmp_v[i].env_->template pool<E>().free(tmp_v[i].slot_);
        tmp_v[i].env_ = env;
        tmp_v[i].slot_ = env->template pool<E>().alloc();
        detail::ad_const(env->layout(), env->template pool<E>(),
                         tmp_v[i].slot_, E(v[i]));
    }
    for (auto& vv : tmp_ov) {
        vv.env_->template pool<E>().free(vv.slot_);
        vv.env_ = env;
        vv.slot_ = env->template pool<E>().alloc();
    }
    std::vector<unsigned> ad_iv, ad_v, ad_ov_slots;
    for (auto& vv : ivecs) ad_iv.push_back(vv.slot_);
    for (auto& vv : tmp_v) ad_v.push_back(vv.slot_);
    for (auto& vv : tmp_ov) ad_ov_slots.push_back(vv.slot_);
    detail::ad_composition(env->layout(), env->template pool<E>(),
                           ad_iv, ad_v, ad_ov_slots);
    // Extract constant term (index 0) from each output slot as Expression
    ovecs.resize(ivecs.size());
    for (unsigned i = 0; i < ivecs.size(); ++i) {
        const E* p = env->template pool<E>().slot(ad_ov_slots[i]);
        ovecs[i] = p[0];
    }
}

#endif // DA_WITH_SYMBOLIC

} // namespace da