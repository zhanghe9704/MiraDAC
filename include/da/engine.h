/**
 * @file engine.h
 * @brief Templated arithmetic kernels in namespace da::detail.
 *
 * @details All kernels take (const Layout&, Pool<T>&, slot_indices...) and
 * operate entirely through raw pointers resolved at the top of each function
 * (per plan A.8 — no env/vector indirection inside inner loops).
 *
 * Algorithms are bit-exact ports of ref/tpsa/src/tpsa.cpp and
 * ref/tpsa/src/tpsa_extend.cc.  The only change is replacing global
 * advec[iv], adveclen[iv], base, prdidx, gnv, gnd, FULL_VEC_LEN, H with
 * pool.slot(iv), pool.len(iv), layout.base(), layout.prdidx(),
 * layout.num_vars(), layout.max_order(), layout.full_len(), layout.H().
 *
 * Explicit instantiation for double is at the bottom of engine.cpp.
 * Expression instantiation (Stage 6) will be added under DA_WITH_SYMBOLIC.
 */

#pragma once

#include "da/layout.h"
#include "da/pool.h"

#include <vector>
#include <complex>
#include <cstddef>

namespace da {
namespace detail {

// ===========================================================================
// Basic copy / reset / const
// ===========================================================================

template<class T>
void ad_copy(Layout& layout, Pool<T>& pool, unsigned isrc, unsigned idst);

template<class T>
void ad_reset(Layout& layout, Pool<T>& pool, unsigned iv);

/// Set slot iv to a constant: all coefficients 0 except [0] = r; length = 1.
template<class T>
void ad_const(Layout& layout, Pool<T>& pool, unsigned iv, T r);

// ===========================================================================
// Arithmetic — in-place (modifying first argument)
// ===========================================================================

/// dst += src   (port of ad_add from tpsa_extend.cc)
template<class T>
void ad_add(Layout& layout, Pool<T>& pool, unsigned idst, unsigned isrc);

/// dst -= src
template<class T>
void ad_sub(Layout& layout, Pool<T>& pool, unsigned idst, unsigned isrc);

/// dst *= c
template<class T>
void ad_mult_const(Layout& layout, Pool<T>& pool, unsigned iv, T c);

/// dst[0] += r  (add a scalar constant to the constant term)
template<class T>
void ad_add_const(Layout& layout, Pool<T>& pool, unsigned iv, T r);

/// iv /= c
template<class T>
void ad_div_c(Layout& layout, Pool<T>& pool, unsigned iv, T c);

// ===========================================================================
// Arithmetic — three-slot (dst must differ from lhs and rhs)
// ===========================================================================

/// dst = lhs * rhs  (port of ad_mult from tpsa_extend.cc)
template<class T>
void ad_mult(Layout& layout, Pool<T>& pool,
             unsigned ilhs, unsigned irhs, unsigned idst);

/// dst = c / iv   (port of ad_c_div from tpsa_extend.cc)
template<class T>
void ad_c_div(Layout& layout, Pool<T>& pool,
              unsigned iv, T c, unsigned ivret);

/// dst = lhs / rhs  (port of ad_div — computes 1/rhs then mult)
template<class T>
void ad_div(Layout& layout, Pool<T>& pool,
            unsigned ilhs, unsigned irhs, unsigned idst);

/// ov = c * iv  (mult_c variant — non-destructive mult by const; port of ad_mult_c)
template<class T>
void ad_mult_c(Layout& layout, Pool<T>& pool,
               unsigned iv, T c, unsigned ov);

// ===========================================================================
// Element access (pok/pek/elem/var)
// ===========================================================================

/// Set coefficient at exponent vector c[] (length n) to value x.
template<class T>
void ad_pok(Layout& layout, Pool<T>& pool,
            unsigned ivec, const int* c, std::size_t n, T x);

/// Get coefficient at exponent vector c[] (length n) into x.
template<class T>
void ad_pek(Layout& layout, Pool<T>& pool,
            unsigned ivec, const int* c, std::size_t n, T& x);

/// Return the idx-th element (1-based) value and exponent array c[].
template<class T>
void ad_elem(Layout& layout, Pool<T>& pool,
             unsigned ivec, unsigned idx, unsigned* c, T& x);

/// Set slot ivec to the ibvec-th base variable with constant term x0.
template<class T>
void ad_var(Layout& layout, Pool<T>& pool,
            unsigned ivec, T x0, unsigned ibvec);

// ===========================================================================
// Nonlinear functions (allocate scratch from pool internally)
// ===========================================================================

template<class T>
void ad_sqrt(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret);

template<class T>
void ad_exp(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret);

template<class T>
void ad_log(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret);

template<class T>
void ad_sin(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret);

template<class T>
void ad_cos(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret);

// ===========================================================================
// Derivative / integrate
// ===========================================================================

/// d(iv)/d(x_{expo}) -> iret
template<class T>
void ad_der(Layout& layout, Pool<T>& pool,
            unsigned iv, unsigned expo, unsigned iret);

/// integral of iv w.r.t. x_{base_id} -> ov
template<class T>
void ad_int(Layout& layout, Pool<T>& pool,
            unsigned iv, unsigned base_id, unsigned ov);

// ===========================================================================
// Composition / substitution
// ===========================================================================

/// TPS composition: ovecs[i] = ivecs[i](v[0],...,v[gnv-1])  (TPS in, TPS out)
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<unsigned>& v,
                    std::vector<unsigned>& ovecs);

/// Composition with double values: ovecs[i] = ivecs[i](v[0],...,v[gnv-1])
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<double>& v,
                    std::vector<double>& ovecs);

/// Composition with complex double values
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<std::complex<double>>& v,
                    std::vector<std::complex<double>>& ovecs);

/// Substitute TPS vectors v into bases base_id[] of ivecs, result in ovecs.
template<class T>
void ad_substitute(Layout& layout, Pool<T>& pool,
                   std::vector<unsigned>& ivecs,
                   std::vector<unsigned>& base_id,
                   std::vector<unsigned>& v,
                   std::vector<unsigned>& ovecs);

/// Single-vector variant: substitute v into base base_id of iv, result in ov.
template<class T>
void ad_substitute(Layout& layout, Pool<T>& pool,
                   unsigned iv, unsigned base_id, unsigned v, unsigned ov);

// ===========================================================================
// Utilities
// ===========================================================================

/// Set coefficients with abs < |eps| to 0; shorten length.
template<class T>
void ad_clean(Layout& layout, Pool<T>& pool, unsigned iv, double eps);

/// Count non-zero coefficients.
template<class T>
int ad_n_element(Layout& layout, Pool<T>& pool, unsigned iv);

/// Return true if all |coeff| < eps (default: double::min).
template<class T>
bool ad_zero_check(Layout& layout, Pool<T>& pool,
                   unsigned iv, double eps = -1.0);

/// Infinity norm (max |coeff|).
template<class T>
double ad_norm(Layout& layout, Pool<T>& pool, unsigned iv);

/// Weighted norm: max |coeff * w^order|.
template<class T>
double ad_weighted_norm(Layout& layout, Pool<T>& pool,
                        unsigned iv, double w);

/// Zero all coefficients; keep length (port of ad_reset_vector).
template<class T>
void ad_reset_vector(Layout& layout, Pool<T>& pool, unsigned iv);

} // namespace detail
} // namespace da
