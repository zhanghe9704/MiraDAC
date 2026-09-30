/**
 * @file engine.cpp
 * @brief Definitions of all da::detail templated kernels + double instantiation.
 *
 * Algorithms are bit-exact ports of:
 *   ref/tpsa/src/tpsa.cpp           (Lingyun Yang)
 *   ref/tpsa/src/tpsa_extend.cc     (He Zhang)
 *
 * Global access pattern in the reference:
 *   advec[iv]       -> pool.slot(iv)
 *   adveclen[iv]    -> pool.len(iv) / pool.set_len(iv,...)
 *   base            -> layout.base()
 *   prdidx          -> layout.prdidx()
 *   gnv             -> layout.num_vars()
 *   gnd             -> layout.max_order()
 *   FULL_VEC_LEN    -> layout.full_len()
 *   H[i][j]         -> layout.H()[i][j]
 *   order_index[k]  -> layout.order_index()[k]
 *
 * Per plan A.8: raw pointers are resolved once at the top of each kernel;
 * inner loops run entirely on those raw pointers.
 */

#include "da/engine.h"

#ifdef DA_WITH_SYMBOLIC
#  include "da/symbolic_ops.h"
#  include <symengine/expression.h>
#  include <symengine/functions.h>
#  include <symengine/eval_double.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>
#include <complex>

namespace da {
namespace detail {

// ---------------------------------------------------------------------------
// Scratch slots of one kernel call, returned to the pool on every exit path.
// A kernel that throws part-way (the pool running out while it takes its next
// temporary) used to leak the temporaries it already held.
// ---------------------------------------------------------------------------
namespace {

// One temporary slot (no heap use: these sit on hot paths).
template<class T>
struct TempSlot {
    Pool<T>& pool;
    unsigned i;
    explicit TempSlot(Pool<T>& p) : pool(p), i(p.alloc()) {}
    ~TempSlot() { pool.free(i); }
    TempSlot(const TempSlot&) = delete;
    TempSlot& operator=(const TempSlot&) = delete;
    operator unsigned() const { return i; }
};

template<class T>
void swap(TempSlot<T>& a, TempSlot<T>& b) noexcept { std::swap(a.i, b.i); }

// A variable number of slots (the power tables).
template<class T>
class ScratchSlots {
public:
    explicit ScratchSlots(Pool<T>& pool) : pool_(pool) {}
    ~ScratchSlots() {
        for (auto it = slots_.rbegin(); it != slots_.rend(); ++it) pool_.free(*it);
    }
    ScratchSlots(const ScratchSlots&) = delete;
    ScratchSlots& operator=(const ScratchSlots&) = delete;

    unsigned assign() {
        slots_.reserve(slots_.size() + 1);  // may throw before a slot is taken
        slots_.push_back(pool_.assign());
        return slots_.back();
    }
    void free(unsigned i) {
        slots_.erase(std::find(slots_.begin(), slots_.end(), i));
        pool_.free(i);
    }

private:
    Pool<T>& pool_;
    std::vector<unsigned> slots_;
};

} // namespace

// ---------------------------------------------------------------------------
// Internal helpers: scalar dispatch for T  (double vs Expression)
// ---------------------------------------------------------------------------
namespace {

// is_near_zero: true if |x| < numeric_limits<double>::min()
// For Expression: uses da::is_zero (symbolic simplify+expand)
template<class T>
bool scalar_is_near_zero(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::abs(x) < std::numeric_limits<double>::min();
    }
#ifdef DA_WITH_SYMBOLIC
    else if constexpr (std::is_same_v<T, SymEngine::Expression>) {
        return da::is_zero(x);
    }
#endif
    else {
        return scalar_is_near_zero(x);
    }
}

// scalar_to_double: extract a double from T
template<class T>
double scalar_to_double(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return x;
    }
#ifdef DA_WITH_SYMBOLIC
    else if constexpr (std::is_same_v<T, SymEngine::Expression>) {
        return SymEngine::eval_double(*x.get_basic());
    }
#endif
    else {
        return static_cast<double>(x);
    }
}

// scalar_sqrt, scalar_exp, scalar_log: scalar function of T returning T
template<class T>
T scalar_sqrt(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::sqrt(x);
    }
#ifdef DA_WITH_SYMBOLIC
    else if constexpr (std::is_same_v<T, SymEngine::Expression>) {
        return SymEngine::Expression(SymEngine::sqrt(x.get_basic()));
    }
#endif
    else {
        return static_cast<T>(std::sqrt(static_cast<double>(x)));
    }
}

template<class T>
T scalar_exp(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::exp(x);
    }
#ifdef DA_WITH_SYMBOLIC
    else if constexpr (std::is_same_v<T, SymEngine::Expression>) {
        return SymEngine::Expression(SymEngine::exp(x.get_basic()));
    }
#endif
    else {
        return static_cast<T>(std::exp(static_cast<double>(x)));
    }
}

template<class T>
T scalar_log(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::log(x);
    }
#ifdef DA_WITH_SYMBOLIC
    else if constexpr (std::is_same_v<T, SymEngine::Expression>) {
        return SymEngine::Expression(SymEngine::log(x.get_basic()));
    }
#endif
    else {
        return static_cast<T>(std::log(static_cast<double>(x)));
    }
}

template<class T>
double scalar_sin(const T& x) {
    return std::sin(scalar_to_double(x));
}

template<class T>
double scalar_cos(const T& x) {
    return std::cos(scalar_to_double(x));
}


// ad_pow_int_pos: internal helper for composition — compute iv^order using
// the power-vector cache power_vv[order], allocating slots as needed.
// power_vv must have size >= gnd+1; power_vv[0] = const-1 slot,
// power_vv[1] = iv slot (NOT owned here), higher indices allocated by first use.
template<class T>
unsigned ad_pow_int_pos_impl(Layout& layout, Pool<T>& pool,
                             unsigned iv,
                             std::vector<unsigned>& power_v,
                             int order, unsigned order_rec)
{
    unsigned res;
    if (pool.len(power_v.at(static_cast<size_t>(order)*order_rec)) > 0) {
        return power_v.at(static_cast<size_t>(order)*order_rec);
    }
    if (order == 0) {
        T* p0 = pool.slot(power_v.at(0));
        p0[0] = T{1};
        pool.set_len(power_v.at(0), 1);
        return power_v.at(0);
    }
    if (order == 1) {
        res = power_v.at(1u * order_rec);
        if (pool.len(res) == 0)
            ad_copy(layout, pool, iv, res);
        return res;
    }
    if (order == 2) {
        res = power_v.at(2u * order_rec);
        if (pool.len(res) == 0)
            ad_mult(layout, pool, iv, iv, res);
        return res;
    }
    if (order & 1) { // odd
        unsigned order_idx = order_rec * 2;
        unsigned vres = power_v.at(order_idx);
        if (pool.len(vres) == 0)
            ad_mult(layout, pool, iv, iv, vres);
        vres = ad_pow_int_pos_impl(layout, pool, vres, power_v, order/2, order_idx);
        res = power_v.at(order_rec + static_cast<size_t>(order/2)*order_idx);
        ad_mult(layout, pool, iv, vres, res);
    } else { // even
        unsigned order_idx = order_rec * 2;
        unsigned vres = power_v.at(order_idx);
        if (pool.len(vres) == 0)
            ad_mult(layout, pool, iv, iv, vres);
        res = ad_pow_int_pos_impl(layout, pool, vres, power_v, order/2, order_idx);
    }
    return res;
}

} // anonymous namespace

// ===========================================================================
// ad_copy
// ===========================================================================
template<class T>
void ad_copy(Layout& layout, Pool<T>& pool, unsigned isrc, unsigned idst)
{
    if (isrc == idst) return;
    unsigned len = pool.len(isrc);
    // Honor a temporarily lowered order (da_change_order). The reference
    // copies FULL_VEC_LEN elements, so a copy made while the order is
    // reduced keeps only the terms up to that order. Without this clamp,
    // "da_change_order(1); t = v; da_restore_order();" does not truncate,
    // and callers such as inv_map see an all-zero nonlinear part.
    const unsigned cap = layout.full_len();
    if (len > cap) len = cap;
    pool.copy_slot(pool.slot(isrc), pool.slot(idst), len);
    pool.set_len(idst, len);
}

// ===========================================================================
// ad_reset  (zero all coefficients, set length = 0)
// ===========================================================================
template<class T>
void ad_reset(Layout& layout, Pool<T>& pool, unsigned iv)
{
    pool.zero_slot(pool.slot(iv));
    pool.set_len(iv, 0);
    (void)layout;
}

// ===========================================================================
// ad_reset_vector  (zero all coefficients; keep length unchanged)
// ===========================================================================
template<class T>
void ad_reset_vector(Layout& layout, Pool<T>& pool, unsigned iv)
{
    // memset full_len elements regardless of current length
    unsigned full_len = layout.full_len();
    T* v = pool.slot(iv);
    if constexpr (std::is_trivially_copyable_v<T>) {
        std::memset(v, 0, full_len * sizeof(T));
    } else {
        for (unsigned k = 0; k < full_len; ++k) v[k] = T{};
    }
}

// ===========================================================================
// ad_const  (set slot to a pure constant r; length = 1)
// ===========================================================================
template<class T>
void ad_const(Layout& /*layout*/, Pool<T>& pool, unsigned iv, T r)
{
    T* v = pool.slot(iv);
    unsigned len = pool.len(iv);
    for (unsigned i = 0; i < len; ++i) v[i] = T{};
    v[0] = r;
    pool.set_len(iv, 1);
}

// ===========================================================================
// ad_add  (dst += src)
// Port of ad_add from tpsa_extend.cc (which replaced the one in tpsa.cpp).
// ===========================================================================
template<class T>
void ad_add(Layout& /*layout*/, Pool<T>& pool, unsigned idst, unsigned isrc)
{
    T*       v    = pool.slot(idst);
    const T* rhsv = pool.slot(isrc);
    unsigned li   = pool.len(idst);
    unsigned lj   = pool.len(isrc);

    if (li < lj) {
        for (unsigned ii = 0; ii < li;  ++ii) v[ii] += rhsv[ii];
        for (unsigned ii = li; ii < lj; ++ii) v[ii]  = rhsv[ii];
        pool.set_len(idst, lj);
    } else {
        for (unsigned ii = 0; ii < lj; ++ii) v[ii] += rhsv[ii];
    }
}

// ===========================================================================
// ad_sub  (dst -= src)
// ===========================================================================
template<class T>
void ad_sub(Layout& /*layout*/, Pool<T>& pool, unsigned idst, unsigned isrc)
{
    T*       v    = pool.slot(idst);
    const T* rhsv = pool.slot(isrc);
    unsigned li   = pool.len(idst);
    unsigned lj   = pool.len(isrc);

    if (li == 0 || lj == 0) return; // reference guard

    if (li < lj) {
        for (unsigned ii = 0; ii < li;  ++ii) v[ii] -= rhsv[ii];
        for (unsigned ii = li; ii < lj; ++ii) v[ii]  = -rhsv[ii];
        pool.set_len(idst, lj);
    } else {
        for (unsigned ii = 0; ii < lj; ++ii) v[ii] -= rhsv[ii];
    }
}

// ===========================================================================
// ad_mult_const  (iv *= c)
// ===========================================================================
template<class T>
void ad_mult_const(Layout& /*layout*/, Pool<T>& pool, unsigned iv, T c)
{
    T*       v = pool.slot(iv);
    unsigned l = pool.len(iv);
    for (unsigned i = 0; i < l; ++i) v[i] *= c;
}

// ===========================================================================
// ad_add_const  (iv[0] += r)
// ===========================================================================
template<class T>
void ad_add_const(Layout& /*layout*/, Pool<T>& pool, unsigned iv, T r)
{
    pool.slot(iv)[0] += r;
}

// ===========================================================================
// ad_div_c  (iv /= c; uses multiply-by-reciprocal like tpsa_extend.cc)
// ===========================================================================
template<class T>
void ad_div_c(Layout& layout, Pool<T>& pool, unsigned iv, T c)
{
    if (scalar_is_near_zero(c)) {
        throw std::invalid_argument("da::detail::ad_div_c: divide by zero");
    }
    const unsigned int* oi = layout.order_index();
    unsigned gnd = layout.max_order();
    T c_inv = T{1} / c;
    T* v = pool.slot(iv);
    unsigned lim = std::min(pool.len(iv), oi[gnd+1]);
    for (unsigned i = 0; i < lim; ++i) v[i] *= c_inv;
}

// ===========================================================================
// ad_mult_c  (ov = c * iv; non-destructive; port of ad_mult_c in tpsa_extend.cc)
// ===========================================================================
template<class T>
void ad_mult_c(Layout& layout, Pool<T>& pool, unsigned iv, T c, unsigned ov)
{
    const unsigned int* oi = layout.order_index();
    unsigned gnd = layout.max_order();
    unsigned len = pool.len(iv);
    unsigned lim = std::min(len, oi[gnd+1]);
    pool.set_len(ov, len);
    const T* src = pool.slot(iv);
    T*       dst = pool.slot(ov);
    for (unsigned i = 0; i < lim; ++i)
        dst[i] = c * src[i];
}

// ===========================================================================
// ad_mult  (dst = lhs * rhs; port of ad_mult from tpsa_extend.cc)
// dst MUST differ from lhs and rhs.
// ===========================================================================
template<class T>
void ad_mult(Layout& layout, Pool<T>& pool,
             unsigned ilhs, unsigned irhs, unsigned idst)
{
    unsigned lhs = ilhs, rhs = irhs, dst = idst;

    const unsigned int* oi     = layout.order_index();
    unsigned int**      pidx   = layout.prdidx();
    unsigned gnd = layout.max_order();
    unsigned full_len = layout.full_len();

    // Zero the destination
    if constexpr (std::is_trivially_copyable_v<T>) {
        std::memset(pool.slot(dst), 0, full_len * sizeof(T));
    } else {
        T* dp = pool.slot(dst);
        for (unsigned k = 0; k < full_len; ++k) dp[k] = T{};
    }

    const T* lv = pool.slot(lhs);
    const T* rv = pool.slot(rhs);
    T*       dv = pool.slot(dst);

    unsigned ll = pool.len(lhs);
    unsigned lr = pool.len(rhs);
    pool.set_len(dst, ll);

    dv[0] = lv[0] * rv[0];

    if (lv[0] != T{}) {
        unsigned lim = std::min(lr, oi[gnd+1]);
        for (unsigned i = 1; i < lim; ++i)
            dv[i] += lv[0] * rv[i];
    }
    if (rv[0] != T{}) {
        unsigned lim = std::min(ll, oi[gnd+1]);
        for (unsigned i = 1; i < lim; ++i)
            dv[i] += lv[i] * rv[0];
    }

    unsigned ord = 1;
    size_t L = std::max(ll, lr);

    unsigned lim_l = std::min(ll, oi[gnd]);
    for (unsigned i = 1; i < lim_l; ++i) {
        if (oi[ord+1] <= i) ++ord;
        unsigned M = oi[gnd-ord+1];
        if (M > lr) M = lr;
        if (lv[i] != T{}) {
            for (unsigned j = 1; j < M; ++j)
                dv[pidx[i][j]] += lv[i] * rv[j];
        }
        if (pidx[i][M-1] >= L) L = pidx[i][M-1] + 1;
    }

    if (L > full_len) L = full_len;
    pool.set_len(dst, static_cast<unsigned>(L));

    // trim trailing zeros
    T*       dv2 = pool.slot(dst);
    unsigned len = pool.len(dst);
    while (len > 1 &&
           scalar_is_near_zero(dv2[len-1]))
        --len;
    pool.set_len(dst, len);
}

// ===========================================================================
// ad_c_div  (ivret = c / iv;  port of ad_c_div from tpsa_extend.cc)
// ===========================================================================
template<class T>
void ad_c_div(Layout& layout, Pool<T>& pool,
              unsigned iv, T c, unsigned ivret)
{
    if (scalar_is_near_zero(c)) {
        ad_reset(layout, pool, ivret);
        return;
    }

    TempSlot<T> ipn(pool);
    TempSlot<T> ip(pool);
    TempSlot<T> itmp(pool);

    unsigned iret = ivret;

    ad_copy(layout, pool, iv, ip);
    ad_copy(layout, pool, iv, ipn);

    T* pn  = pool.slot(ipn);
    T* p   = pool.slot(ip);
    T* ret = pool.slot(iret);
    T* v   = pool.slot(iv);

    T x0 = v[0];
    T inv_x0 = T{1} / x0;
    inv_x0 = -inv_x0;

    unsigned lp = pool.len(ip);
    for (unsigned i = 1; i < lp; ++i) p[i] *= inv_x0;
    pn[0] = p[0] = T{};
    inv_x0 = -inv_x0; // restore positive

    // Copy p into pn (they should already match after ad_copy, just update [0])
    if constexpr (std::is_trivially_copyable_v<T>) {
        std::memcpy(pn, p, lp * sizeof(T));
    } else {
        for (unsigned k = 0; k < lp; ++k) pn[k] = p[k];
    }

    ret[0] = T{1};
    pool.set_len(iret, 1);

    ad_add(layout, pool, iret, ip);

    unsigned gnd = layout.max_order();
    unsigned full_len = layout.full_len();
    for (unsigned nd = 2; nd < gnd+1; ++nd) {
        ad_mult(layout, pool, ipn, ip, itmp);
        ad_copy(layout, pool, itmp, ipn);
        pn = pool.slot(ipn);
        unsigned lm = pool.len(ipn);
        ret = pool.slot(iret);
        for (unsigned i = 0; i < lm; ++i) ret[i] += pn[i];
    }
    pool.set_len(iret, full_len);

    T ret_coef = inv_x0 * c;
    ret = pool.slot(iret);
    unsigned lr = pool.len(iret);
    for (unsigned i = 0; i < lr; ++i) ret[i] *= ret_coef;

}

// ===========================================================================
// ad_div  (dst = lhs / rhs)
// ===========================================================================
template<class T>
void ad_div(Layout& layout, Pool<T>& pool,
            unsigned ilhs, unsigned irhs, unsigned idst)
{
    T c = T{1};
    TempSlot<T> itmp(pool);
    ad_c_div(layout, pool, irhs, c, itmp);
    ad_mult(layout, pool, ilhs, itmp, idst);
}

// ===========================================================================
// ad_pok  (set coefficient at exponent vector)
// ===========================================================================
template<class T>
void ad_pok(Layout& layout, Pool<T>& pool,
            unsigned ivec, const int* c, std::size_t n, T x)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const auto& H = layout.H();

    std::size_t N = (n > gnv) ? gnv : n;
    unsigned d = 0;
    std::vector<unsigned> cef(gnv, 0);
    for (std::size_t i = 0; i < N; ++i) {
        cef[i] = static_cast<unsigned>(c[i]);
        d += cef[i];
    }

    if (d > gnd) return;

    std::size_t k = 0;
    for (unsigned i = 0; i < gnv; ++i) {
        unsigned bv_i = d;
        d -= cef[i];
        k += H[gnv-i][bv_i];
    }

    T* v = pool.slot(ivec);
    v[k] = x;
    if (k + 1 > pool.len(ivec))
        pool.set_len(ivec, static_cast<unsigned>(k + 1));
}

// ===========================================================================
// ad_pek  (get coefficient at exponent vector)
// ===========================================================================
template<class T>
void ad_pek(Layout& layout, Pool<T>& pool,
            unsigned ivec, const int* c, std::size_t n, T& x)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const auto& H = layout.H();

    std::size_t N = (n > gnv) ? gnv : n;
    unsigned d = 0;
    std::vector<unsigned> cef(gnv, 0);
    for (std::size_t i = 0; i < N; ++i) {
        cef[i] = static_cast<unsigned>(c[i]);
        d += cef[i];
    }

    if (d > gnd) { x = T{}; return; }

    std::size_t k = 0;
    for (unsigned i = 0; i < gnv; ++i) {
        unsigned bv_i = d;
        d -= cef[i];
        k += H[gnv-i][bv_i];
    }

    if (k > pool.len(ivec)) x = T{};
    else x = pool.slot(ivec)[k];
}

// ===========================================================================
// ad_elem  (idx-th element, 1-based)
// ===========================================================================
template<class T>
void ad_elem(Layout& layout, Pool<T>& pool,
             unsigned ivec, unsigned idx, unsigned* c, T& x)
{
    unsigned gnv = layout.num_vars();
    const unsigned int* bptr = layout.base();

    for (unsigned i = 0; i < gnv; ++i) c[i] = 0;
    if (idx > pool.len(ivec) || idx < 1) { x = T{}; return; }

    x = pool.slot(ivec)[idx - 1];
    const unsigned int* p = bptr + static_cast<std::size_t>(gnv) * (idx - 1);
    for (unsigned j = 0; j < gnv - 1; ++j) {
        c[j] = *p - *(p+1);
        ++p;
    }
    c[gnv-1] = *p;
}

// ===========================================================================
// ad_var  (set slot to a base variable)
// ===========================================================================
template<class T>
void ad_var(Layout& layout, Pool<T>& pool,
            unsigned ivec, T x0, unsigned ibvec)
{
    unsigned gnv = layout.num_vars();
    unsigned full_len = layout.full_len();
    T* v = pool.slot(ivec);

    for (unsigned k = 0; k < full_len; ++k) v[k] = T{};
    v[0] = x0;

    if (ibvec < gnv) {
        pool.set_len(ivec, ibvec + 2);
        v[ibvec + 1] = T{1};
    } else {
        pool.set_len(ivec, 1);
    }
}

// ===========================================================================
// ad_sqrt  (port of ad_sqrt from tpsa.cpp)
// ===========================================================================
template<class T>
void ad_sqrt(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret)
{
    T x = pool.slot(iv)[0];
    unsigned gnd = layout.max_order();

    TempSlot<T> itmp(pool);
    TempSlot<T> ip(pool);
    TempSlot<T> ipn(pool);

    ad_copy(layout, pool, iv, ip);
    ad_div_c(layout, pool, ip, x);
    pool.slot(ip)[0] = T{};

    ad_reset(layout, pool, iret);
    pool.slot(iret)[0] = T{1};
    pool.set_len(iret, 1);

    ad_copy(layout, pool, ip, itmp);
    ad_copy(layout, pool, ip, ipn);

    double c = 0.5;
    for (unsigned i = 1; i < gnd + 1; ++i) {
        ad_mult_const(layout, pool, itmp, T{c});
        ad_add(layout, pool, iret, itmp);
        c = c * (1.0 - 2.0*i) / 2.0 / (i + 1.0);
        ad_mult(layout, pool, ip, ipn, itmp);
        ad_copy(layout, pool, itmp, ipn);
    }

    T sx = scalar_sqrt(x);
    ad_mult_const(layout, pool, iret, sx);

}

// ===========================================================================
// ad_exp  (port of ad_exp from tpsa.cpp)
// ===========================================================================
template<class T>
void ad_exp(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret)
{
    T ex = scalar_exp(pool.slot(iv)[0]);
    unsigned gnd = layout.max_order();

    TempSlot<T> itmp(pool);
    TempSlot<T> ip(pool);
    TempSlot<T> ipn(pool);

    ad_copy(layout, pool, iv, ip);
    pool.slot(ip)[0] = T{};

    ad_reset(layout, pool, iret);
    pool.slot(iret)[0] = T{1};
    pool.set_len(iret, 1);

    ad_copy(layout, pool, ip, itmp);
    ad_copy(layout, pool, ip, ipn);

    double c = 1.0;
    for (unsigned i = 1; i < gnd + 1; ++i) {
        c = c * static_cast<double>(i);
        ad_div_c(layout, pool, itmp, T{c});
        ad_add(layout, pool, iret, itmp);
        ad_mult(layout, pool, ip, ipn, itmp);
        ad_copy(layout, pool, itmp, ipn);
    }

    ad_mult_const(layout, pool, iret, ex);

}

// ===========================================================================
// ad_log  (port of ad_log from tpsa.cpp)
// ===========================================================================
template<class T>
void ad_log(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret)
{
    T x0   = pool.slot(iv)[0];
    T logx0 = scalar_log(x0);
    unsigned gnd = layout.max_order();

    TempSlot<T> itmp(pool);
    TempSlot<T> ip(pool);
    TempSlot<T> ipn(pool);

    ad_copy(layout, pool, iv, ip);
    ad_div_c(layout, pool, ip, x0);
    pool.slot(ip)[0] = T{};

    ad_reset(layout, pool, iret);
    pool.slot(iret)[0] = logx0;
    pool.set_len(iret, 1);

    ad_copy(layout, pool, ip, itmp);
    ad_copy(layout, pool, ip, ipn);

    for (unsigned i = 1; i < gnd + 1; ++i) {
        double c = (i % 2 == 0) ? -1.0*i : 1.0*i;
        ad_div_c(layout, pool, itmp, T{c});
        ad_add(layout, pool, iret, itmp);
        ad_mult(layout, pool, ip, ipn, itmp);
        ad_copy(layout, pool, itmp, ipn);
    }

}

// ===========================================================================
// ad_sin  (port of ad_sin from tpsa.cpp)
// ===========================================================================
template<class T>
void ad_sin(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret)
{
    unsigned gnv      = layout.num_vars();
    unsigned gnd      = layout.max_order();
    unsigned full_len = layout.full_len();
    (void)gnv;

    T* v0 = pool.slot(iv);
    double s = scalar_sin(v0[0]);
    double c = scalar_cos(v0[0]);

    TempSlot<T> ipnev(pool);
    TempSlot<T> ipnod(pool);
    TempSlot<T> ip(pool);

    ad_copy(layout, pool, iv, iret);
    ad_copy(layout, pool, iv, ipnev);
    ad_copy(layout, pool, iv, ipnod);
    ad_copy(layout, pool, iv, ip);

    T* ret  = pool.slot(iret);
    T* pnod = pool.slot(ipnod);
    T* pnev = pool.slot(ipnev);
    T* p    = pool.slot(ip);

    unsigned plen = pool.len(ip);

    pnev[0] = pnod[0] = p[0] = T{};
    ret[0] = static_cast<T>(s);

    for (unsigned i = 1; i < plen; ++i)
        ret[i] *= static_cast<T>(c);

    for (unsigned k = 2; k < gnd + 1; ++k) {
        ad_mult(layout, pool, ip, ipnod, ipnev);
        pnev = pool.slot(ipnev);
        unsigned lev = pool.len(ipnev);
        ret  = pool.slot(iret);
        for (unsigned i = 0; i < lev; ++i) {
            pnev[i] /= static_cast<T>(static_cast<double>(k));
            switch (k % 4) {
            case 0: ret[i] += static_cast<T>(s) * pnev[i]; break;
            case 1: ret[i] += static_cast<T>(c) * pnev[i]; break;
            case 2: ret[i] -= static_cast<T>(s) * pnev[i]; break;
            case 3: ret[i] -= static_cast<T>(c) * pnev[i]; break;
            }
        }
        swap(ipnev, ipnod);
        pnev = pool.slot(ipnev);
        pnod = pool.slot(ipnod);
    }

    pool.set_len(iret, full_len);

}

// ===========================================================================
// ad_cos  (port of ad_cos from tpsa.cpp)
// ===========================================================================
template<class T>
void ad_cos(Layout& layout, Pool<T>& pool, unsigned iv, unsigned iret)
{
    unsigned gnd      = layout.max_order();
    unsigned full_len = layout.full_len();

    T* v0 = pool.slot(iv);
    double s = scalar_sin(v0[0]);
    double c = scalar_cos(v0[0]);

    TempSlot<T> ipnev(pool);
    TempSlot<T> ipnod(pool);
    TempSlot<T> ip(pool);

    ad_copy(layout, pool, iv, iret);
    ad_copy(layout, pool, iv, ipnev);
    ad_copy(layout, pool, iv, ipnod);
    ad_copy(layout, pool, iv, ip);

    T* ret  = pool.slot(iret);
    T* pnev = pool.slot(ipnev);
    T* pnod = pool.slot(ipnod);
    T* p    = pool.slot(ip);

    unsigned plen = pool.len(ip);

    pnev[0] = pnod[0] = p[0] = T{};
    ret[0] = static_cast<T>(c);

    for (unsigned i = 1; i < plen; ++i)
        ret[i] *= static_cast<T>(-s);

    for (unsigned k = 2; k < gnd + 1; ++k) {
        ad_mult(layout, pool, ip, ipnod, ipnev);
        pnev = pool.slot(ipnev);
        unsigned lev = pool.len(ipnev);
        ret  = pool.slot(iret);
        for (unsigned i = 0; i < lev; ++i) {
            pnev[i] /= static_cast<T>(static_cast<double>(k));
            switch (k % 4) {
            case 0: ret[i] += static_cast<T>(c) * pnev[i]; break;
            case 1: ret[i] -= static_cast<T>(s) * pnev[i]; break;
            case 2: ret[i] -= static_cast<T>(c) * pnev[i]; break;
            case 3: ret[i] += static_cast<T>(s) * pnev[i]; break;
            }
        }
        swap(ipnev, ipnod);
        pnev = pool.slot(ipnev);
        pnod = pool.slot(ipnod);
    }

    pool.set_len(iret, full_len);

}

// ===========================================================================
// ad_der  (derivative;  port of ad_der from tpsa_extend.cc)
// ===========================================================================
template<class T>
void ad_der(Layout& layout, Pool<T>& pool,
            unsigned iv, unsigned expo, unsigned iret)
{
    unsigned gnv = layout.num_vars();
    const unsigned int* bptr = layout.base();
    const auto& H = layout.H();

    ad_reset(layout, pool, iret);
    pool.slot(iret)[0] = T{};
    pool.set_len(iret, 1);

    const unsigned int* p = bptr;
    unsigned l = pool.len(iv);
    const T* vsrc = pool.slot(iv);
    T* vret = pool.slot(iret);

    std::vector<unsigned> cef(gnv);
    std::vector<unsigned> bv(gnv);

    for (unsigned i = 0; i < l; ++i) {
        unsigned d = 0;
        for (unsigned j = 0; j < gnv - 1; ++j) {
            cef[j] = *p - *(p+1);
            ++p;
            d += cef[j];
        }
        cef[gnv-1] = *p;
        d += *p;
        ++p;

        if (cef[expo] == 0) continue;

        unsigned jexp = cef[expo];
        cef[expo] -= 1;
        --d;

        std::size_t k = 0;
        unsigned dtmp = d;
        for (unsigned j = 0; j < gnv; ++j) {
            bv[j] = dtmp;
            dtmp -= cef[j];
            k += H[gnv-j][bv[j]];
        }

        vret[k] = vsrc[i] * static_cast<T>(static_cast<double>(jexp));
        if (k >= pool.len(iret))
            pool.set_len(iret, static_cast<unsigned>(k + 1));
    }
}

// ===========================================================================
// ad_int  (integrate; port of ad_int from tpsa_extend.cc)
// ===========================================================================
template<class T>
void ad_int(Layout& layout, Pool<T>& pool,
            unsigned iv, unsigned base_id, unsigned ov)
{
    unsigned gnv = layout.num_vars();
    const unsigned int* bptr = layout.base();

    // multiply by the base variable, then divide each coefficient by its order
    TempSlot<T> vtemp(pool);
    T x0 = T{};
    ad_var(layout, pool, vtemp, x0, base_id);
    ad_mult(layout, pool, iv, vtemp, ov);

    // Now divide each coefficient by the exponent of base_id in that monomial
    const unsigned int* p = bptr;
    unsigned lv = pool.len(ov);
    T* vov = pool.slot(ov);

    for (unsigned i = 0; i < lv; ++i) {
        if (scalar_is_near_zero(vov[i])) {
            p += gnv;
            continue;
        }
        // decode exponent of base_id
        std::vector<unsigned> c(gnv);
        const unsigned int* pp = bptr + static_cast<std::size_t>(gnv) * i;
        for (unsigned j = 0; j < gnv - 1; ++j) {
            c[j] = *pp - *(pp+1);
            ++pp;
        }
        c[gnv-1] = *pp;
        // Use integer cast (not double) to preserve exact rational arithmetic
        // when T=SymEngine::Expression (avoids floating-point coefficients).
        vov[i] /= static_cast<T>(static_cast<int>(c[base_id]));
    }
    (void)p;

}

// ===========================================================================
// ad_composition (TPS in, TPS out) — port of ad_composition(ivecs, v (TVEC), ovecs)
// from tpsa_extend.cc
// ===========================================================================
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<unsigned>& v,
                    std::vector<unsigned>& ovecs)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    unsigned full_len = layout.full_len();
    const unsigned int* bptr = layout.base();
    const unsigned int* oi   = layout.order_index();
    unsigned int** pidx      = layout.prdidx();

    assert(gnv == v.size());
    assert(ivecs.size() == ovecs.size());

    // Allocate power table: power_vv[var][order] = slot index
    std::vector<std::vector<unsigned>> power_vv(gnv, std::vector<unsigned>(gnd+1));
    ScratchSlots<T> powers(pool);
    for (unsigned var = 0; var < gnv; ++var) {
        for (unsigned ord = 0; ord < gnd+1; ++ord) {
            unsigned idx = powers.assign();
            power_vv[var][ord] = idx;
        }
        // power[0] = constant 1
        T* p0 = pool.slot(power_vv[var][0]);
        p0[0] = T{1};
        pool.set_len(power_vv[var][0], 1);
    }

    // power[1] points to the input variable (borrow slot, free it at end)
    for (unsigned i = 0; i < gnv; ++i) {
        powers.free(power_vv[i][1]);
        power_vv[i][1] = v[i];
    }

    TempSlot<T> tmp(pool);
    TempSlot<T> product(pool);

    // Find max length of input vectors
    unsigned veclen_max = 0;
    for (auto iv : ivecs)
        if (pool.len(iv) > veclen_max) veclen_max = pool.len(iv);

    // Reset output vectors and copy constant element
    for (auto ov : ovecs) {
        ad_reset(layout, pool, ov);
        pool.set_len(ov, 1);
    }
    for (unsigned i = 0; i < ivecs.size(); ++i)
        pool.slot(ovecs[i])[0] = pool.slot(ivecs[i])[0];

    std::vector<unsigned> c(gnv);
    const unsigned int* p = bptr;
    // Read index 0 (constant monomial) c values
    for (unsigned j = 0; j < gnv-1; ++j) {
        c[j] = *p - *(p+1); ++p;
    }
    c[gnv-1] = *p++;

    unsigned vec_size = static_cast<unsigned>(ivecs.size());

    for (unsigned i = 1; i < veclen_max; ++i) {
        bool c_flag = true;
        bool product_flag = true;
        unsigned zero_coef = 0;

        for (unsigned iv = 0; iv < vec_size; ++iv) {
            if (i >= pool.len(ivecs[iv])) { ++zero_coef; continue; }
            if (scalar_is_near_zero(pool.slot(ivecs[iv])[i])) { ++zero_coef; continue; }

            if (c_flag) {
                for (unsigned j = 0; j < gnv-1; ++j) {
                    c[j] = *p - *(p+1); ++p;
                }
                c[gnv-1] = *p++;

                for (unsigned id = 0; id < gnv; ++id) {
                    if (c[id] > 0) {
                        ad_pow_int_pos_impl(layout, pool, v[id], power_vv[id],
                                            static_cast<int>(c[id]), 1);
                    }
                }
                c_flag = false;
            }

            T coef = pool.slot(ivecs[iv])[i];
            if (product_flag) {
                ad_reset_vector(layout, pool, product);
                pool.slot(product)[0] = T{1};
                pool.set_len(product, 1);
                for (unsigned id = 0; id < gnv; ++id) {
                    if (c[id] > 0) {
                        ad_mult(layout, pool, product,
                                power_vv[id][c[id]], tmp);
                        swap(product, tmp);
                    }
                }
                product_flag = false;
            }

            ad_mult_c(layout, pool, product, coef, tmp);
            T* ov_ptr = pool.slot(ovecs[iv]);
            const T* tp = pool.slot(tmp);
            unsigned lt = pool.len(tmp);
            for (unsigned idx = 0; idx < lt; ++idx)
                ov_ptr[idx] += tp[idx];
        }

        if (zero_coef == vec_size) p += gnv;
    }

    // Compute final lengths
    for (auto ov : ovecs) {
        unsigned len = 1;
        T* ov_ptr = pool.slot(ov);
        for (int i = static_cast<int>(oi[gnd+1]) - 1; i >= 0; --i) {
            if (!scalar_is_near_zero(ov_ptr[i])) {
                len = static_cast<unsigned>(i) + 1;
                break;
            }
        }
        pool.set_len(ov, len);
    }


    // power_vv[var][1] == v[var] is borrowed; `powers` frees the rest.

    (void)full_len; (void)pidx;
}

// ===========================================================================
// ad_composition (TPS in, double values out)
// ===========================================================================
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<double>& v,
                    std::vector<double>& ovecs)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const unsigned int* bptr = layout.base();

    assert(gnv == v.size());
    assert(ivecs.size() == ovecs.size());

    // Precompute powers: power_vv[var][ord] = v[var]^ord
    std::vector<std::vector<double>> power_vv(gnv, std::vector<double>(gnd+1));
    for (auto& pv : power_vv) pv[0] = 1.0;
    for (unsigned idx = 0; idx < gnv; ++idx) {
        double val = v[idx];
        for (unsigned i = 1; i < gnd+1; ++i)
            power_vv[idx][i] = val * power_vv[idx][i-1];
    }

    unsigned veclen_max = 0;
    for (auto iv : ivecs)
        if (pool.len(iv) > veclen_max) veclen_max = pool.len(iv);

    for (unsigned i = 0; i < ovecs.size(); ++i)
        ovecs[i] = static_cast<double>(pool.slot(ivecs[i])[0]);

    std::vector<unsigned> c(gnv);
    const unsigned int* p = bptr;
    for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
    c[gnv-1] = *p++;

    unsigned vec_size = static_cast<unsigned>(ivecs.size());
    for (unsigned i = 1; i < veclen_max; ++i) {
        bool product_flag = true;
        bool c_flag = true;
        double product = 1.0;
        unsigned zero_coef = 0;

        for (unsigned iv = 0; iv < vec_size; ++iv) {
            if (i >= pool.len(ivecs[iv])) { ++zero_coef; continue; }
            double coef = static_cast<double>(pool.slot(ivecs[iv])[i]);
            if (std::abs(coef) < std::numeric_limits<double>::min()) {
                ++zero_coef; continue;
            }

            if (c_flag) {
                for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
                c[gnv-1] = *p++;
                c_flag = false;
            }
            if (product_flag) {
                product = 1.0;
                for (unsigned id = 0; id < gnv; ++id)
                    product *= power_vv[id][c[id]];
                product_flag = false;
            }
            ovecs[iv] += product * coef;
        }
        if (zero_coef == vec_size) p += gnv;
    }
}

// ===========================================================================
// ad_composition (TPS in, complex<double> values out)
// ===========================================================================
template<class T>
void ad_composition(Layout& layout, Pool<T>& pool,
                    std::vector<unsigned>& ivecs,
                    std::vector<std::complex<double>>& v,
                    std::vector<std::complex<double>>& ovecs)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const unsigned int* bptr = layout.base();

    assert(gnv == v.size());
    assert(ivecs.size() == ovecs.size());

    std::vector<std::vector<std::complex<double>>> power_vv(gnv,
        std::vector<std::complex<double>>(gnd+1));
    for (auto& pv : power_vv) pv[0] = 1.0;
    for (unsigned idx = 0; idx < gnv; ++idx) {
        auto val = v[idx];
        for (unsigned i = 1; i < gnd+1; ++i)
            power_vv[idx][i] = val * power_vv[idx][i-1];
    }

    unsigned veclen_max = 0;
    for (auto iv : ivecs)
        if (pool.len(iv) > veclen_max) veclen_max = pool.len(iv);

    for (unsigned i = 0; i < ovecs.size(); ++i)
        ovecs[i] = static_cast<double>(pool.slot(ivecs[i])[0]);

    std::vector<unsigned> c(gnv);
    const unsigned int* p = bptr;
    for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
    c[gnv-1] = *p++;

    unsigned vec_size = static_cast<unsigned>(ivecs.size());
    for (unsigned i = 1; i < veclen_max; ++i) {
        bool product_flag = true;
        bool c_flag = true;
        std::complex<double> product = 1.0;
        unsigned zero_coef = 0;

        for (unsigned iv = 0; iv < vec_size; ++iv) {
            if (i >= pool.len(ivecs[iv])) { ++zero_coef; continue; }
            double coef = static_cast<double>(pool.slot(ivecs[iv])[i]);
            if (std::abs(coef) < std::numeric_limits<double>::min()) {
                ++zero_coef; continue;
            }

            if (c_flag) {
                for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
                c[gnv-1] = *p++;
                c_flag = false;
            }
            if (product_flag) {
                product = 1.0;
                for (unsigned id = 0; id < gnv; ++id)
                    product *= power_vv[id][c[id]];
                product_flag = false;
            }
            ovecs[iv] += product * coef;
        }
        if (zero_coef == vec_size) p += gnv;
    }
}

// ===========================================================================
// ad_substitute (vector into single base of single TPS)
// ===========================================================================
template<class T>
void ad_substitute(Layout& layout, Pool<T>& pool,
                   unsigned iv, unsigned base_id, unsigned v_slot, unsigned ov)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const unsigned int* bptr = layout.base();
    const unsigned int* oi   = layout.order_index();
    unsigned int** pidx      = layout.prdidx();
    const auto& H            = layout.H();

    // Build power table for v_slot
    std::vector<unsigned> power_v(gnd+1);
    ScratchSlots<T> powers(pool);
    for (auto& idx : power_v) {
        unsigned si = powers.assign();
        idx = si;
    }
    T* p0 = pool.slot(power_v[0]);
    p0[0] = T{1};
    pool.set_len(power_v[0], 1);
    ad_copy(layout, pool, v_slot, power_v[1]);

    ad_reset(layout, pool, ov);

    unsigned l = pool.len(iv);
    const T*          vsrc = pool.slot(iv);
    const unsigned int* p  = bptr;

    std::vector<unsigned> c(gnv);
    std::vector<unsigned> bv(gnv);

    for (unsigned i = 0; i < l; ++i) {
        if (scalar_is_near_zero(vsrc[i])) {
            p += gnv;
            continue;
        }
        for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
        c[gnv-1] = *p++;

        if (c[base_id] > 0) {
            unsigned tmp_slot = ad_pow_int_pos_impl(layout, pool, v_slot, power_v,
                                                    static_cast<int>(c[base_id]), 1);
            T coef = vsrc[i];
            c[base_id] = 0;

            unsigned d = 0;
            for (unsigned k = 0; k < gnv; ++k) d += c[k];
            unsigned order = gnd - d;
            unsigned idx_limit = oi[order+1];

            unsigned k = 0;
            unsigned dtmp = d;
            for (unsigned k2 = 0; k2 < gnv; ++k2) {
                bv[k2] = dtmp;
                dtmp -= c[k2];
                k += H[gnv-k2][bv[k2]];
            }

            T* ov_ptr = pool.slot(ov);
            const T* tp = pool.slot(tmp_slot);
            unsigned lt = pool.len(tmp_slot);

            if (k > 0) {
                if (!scalar_is_near_zero(tp[0]))
                    ov_ptr[k] += coef * tp[0];
                for (unsigned idx = 1; idx < lt && idx < idx_limit; ++idx)
                    ov_ptr[pidx[k][idx]] += coef * tp[idx];
            } else {
                for (unsigned idx = 0; idx < lt; ++idx) {
                    if (!scalar_is_near_zero(tp[idx]))
                        ov_ptr[idx] += coef * tp[idx];
                }
            }
        } else {
            pool.slot(ov)[i] += vsrc[i];
        }
    }

    // Compute final length
    unsigned len = 1;
    T* ov_ptr = pool.slot(ov);
    for (int i = static_cast<int>(oi[gnd+1]) - 1; i >= 0; --i) {
        if (!scalar_is_near_zero(ov_ptr[i])) {
            len = static_cast<unsigned>(i) + 1;
            break;
        }
    }
    pool.set_len(ov, len);

}

// ===========================================================================
// ad_substitute (multi-base group version)
// ===========================================================================
template<class T>
void ad_substitute(Layout& layout, Pool<T>& pool,
                   std::vector<unsigned>& ivecs,
                   std::vector<unsigned>& base_id,
                   std::vector<unsigned>& v,
                   std::vector<unsigned>& ovecs)
{
    assert(base_id.size() == v.size());
    assert(ivecs.size() == ovecs.size());

    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const unsigned int* bptr = layout.base();
    const unsigned int* oi   = layout.order_index();
    unsigned int** pidx      = layout.prdidx();
    const auto& H            = layout.H();
    unsigned nv = static_cast<unsigned>(v.size());

    // Build power tables
    std::vector<std::vector<unsigned>> power_vv(nv, std::vector<unsigned>(gnd+1));
    ScratchSlots<T> powers(pool);
    for (unsigned k = 0; k < nv; ++k) {
        for (auto& idx : power_vv[k]) idx = powers.assign();
        T* p0 = pool.slot(power_vv[k][0]);
        p0[0] = T{1};
        pool.set_len(power_vv[k][0], 1);
    }
    for (unsigned i = 0; i < nv; ++i) {
        powers.free(power_vv[i][1]);
        power_vv[i][1] = v[i];
    }

    for (auto ov : ovecs) ad_reset(layout, pool, ov);

    TempSlot<T> tmp_slot(pool);
    TempSlot<T> product(pool);

    std::vector<unsigned> c(gnv);
    std::vector<unsigned> bv(gnv);
    std::vector<unsigned> rc(nv, 0);

    unsigned veclen_max = 0;
    for (auto iv : ivecs) if (pool.len(iv) > veclen_max) veclen_max = pool.len(iv);

    const unsigned int* p = bptr + gnv; // skip index 0

    for (unsigned i = 1; i < veclen_max; ++i) {
        bool c_flag = true;
        bool k_flag = true;
        bool sub_flag = false;
        bool product_flag = true;
        unsigned i_count = 0;
        std::fill(rc.begin(), rc.end(), 0u);
        unsigned k = 0;
        unsigned order = 0, idx_limit = 0;

        for (unsigned ivv = 0; ivv < ivecs.size(); ++ivv) {
            if (i >= pool.len(ivecs[ivv])) { ++i_count; continue; }
            if (scalar_is_near_zero(pool.slot(ivecs[ivv])[i])) { ++i_count; continue; }

            if (c_flag) {
                for (unsigned j = 0; j < gnv-1; ++j) { c[j] = *p - *(p+1); ++p; }
                c[gnv-1] = *p++;

                for (unsigned id = 0; id < nv; ++id) {
                    if (c[base_id[id]] > 0) {
                        ad_pow_int_pos_impl(layout, pool, v[id], power_vv[id],
                                            static_cast<int>(c[base_id[id]]), 1);
                        rc[id] = c[base_id[id]];
                        c[base_id[id]] = 0;
                        sub_flag = true;
                    }
                }
                c_flag = false;
            }

            if (sub_flag) {
                if (k_flag) {
                    unsigned d = 0;
                    for (unsigned jj = 0; jj < gnv; ++jj) d += c[jj];
                    order = gnd - d;
                    idx_limit = oi[order+1];
                    unsigned dtmp = d;
                    k = 0;
                    for (unsigned jj = 0; jj < gnv; ++jj) {
                        bv[jj] = dtmp;
                        dtmp -= c[jj];
                        k += H[gnv-jj][bv[jj]];
                    }
                    k_flag = false;
                }

                T coef = pool.slot(ivecs[ivv])[i];
                if (k > 0) {
                    if (product_flag) {
                        ad_reset_vector(layout, pool, product);
                        pool.slot(product)[0] = T{1};
                        pool.set_len(product, 1);
                        layout.change_order(order); // temporarily lower order
                        for (unsigned id = 0; id < nv; ++id) {
                            if (rc[id] > 0) {
                                ad_mult(layout, pool, product, power_vv[id][rc[id]], tmp_slot);
                                swap(product, tmp_slot);
                            }
                        }
                        layout.change_order(gnd);   // back to the caller's order
                        product_flag = false;
                    }
                    ad_mult_c(layout, pool, product, coef, tmp_slot);
                    T* ov_ptr = pool.slot(ovecs[ivv]);
                    const T* tp = pool.slot(tmp_slot);
                    if (!scalar_is_near_zero(tp[0]))
                        ov_ptr[k] += tp[0];
                    unsigned lt = pool.len(tmp_slot);
                    for (unsigned idx = 1; idx < lt && idx < idx_limit; ++idx)
                        ov_ptr[pidx[k][idx]] += tp[idx];
                } else {
                    if (product_flag) {
                        ad_reset_vector(layout, pool, product);
                        pool.slot(product)[0] = T{1};
                        pool.set_len(product, 1);
                        for (unsigned id = 0; id < nv; ++id) {
                            if (rc[id] > 0) {
                                ad_mult(layout, pool, product, power_vv[id][rc[id]], tmp_slot);
                                swap(product, tmp_slot);
                            }
                        }
                        product_flag = false;
                    }
                    ad_mult_c(layout, pool, product, coef, tmp_slot);
                    T* ov_ptr = pool.slot(ovecs[ivv]);
                    const T* tp = pool.slot(tmp_slot);
                    unsigned lt = pool.len(tmp_slot);
                    for (unsigned idx = 0; idx < lt; ++idx)
                        ov_ptr[idx] += tp[idx];
                }
            } else {
                pool.slot(ovecs[ivv])[i] += pool.slot(ivecs[ivv])[i];
            }
        }
        if (static_cast<unsigned>(i_count) == static_cast<unsigned>(ivecs.size())) p += gnv;
    }

    // Add constant terms
    for (unsigned i = 0; i < ivecs.size(); ++i)
        pool.slot(ovecs[i])[0] += pool.slot(ivecs[i])[0];

    for (auto ov : ovecs) {
        unsigned len = 1;
        T* ov_ptr = pool.slot(ov);
        for (int i = static_cast<int>(oi[gnd+1]) - 1; i >= 0; --i) {
            if (!scalar_is_near_zero(ov_ptr[i])) {
                len = static_cast<unsigned>(i) + 1;
                break;
            }
        }
        pool.set_len(ov, len);
    }

    // power_vv[k][1] == v[k] is borrowed; `powers` frees the rest.
}

// ===========================================================================
// ad_clean
// ===========================================================================
template<class T>
void ad_clean(Layout& /*layout*/, Pool<T>& pool, unsigned iv, double eps)
{
    T* p = pool.slot(iv);
    unsigned l = pool.len(iv);
    unsigned N = 0;
    double abseps = std::abs(eps);
    for (unsigned i = 0; i < l; ++i) {
#ifdef DA_WITH_SYMBOLIC
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            if (da::is_zero(p[i]))
                p[i] = T{};
            else
                N = i;
        } else
#endif
        {
            if (std::abs(p[i]) < abseps)
                p[i] = T{};
            else
                N = i;
        }
    }
    if (l > N + 1) pool.set_len(iv, N + 1);
}

// ===========================================================================
// ad_n_element
// ===========================================================================
template<class T>
int ad_n_element(Layout& /*layout*/, Pool<T>& pool, unsigned iv)
{
    int n = 0;
    const T* p = pool.slot(iv);
    unsigned l = pool.len(iv);
    for (unsigned i = 0; i < l; ++i)
        if (!scalar_is_near_zero(p[i]))
            ++n;
    return n;
}

// ===========================================================================
// ad_zero_check
// ===========================================================================
template<class T>
bool ad_zero_check(Layout& /*layout*/, Pool<T>& pool,
                   unsigned iv, double eps)
{
    double zero = std::numeric_limits<double>::min();
    if (eps > 0) zero = eps;
    const T* p = pool.slot(iv);
    unsigned l = pool.len(iv);
    for (unsigned i = 0; i < l; ++i)
#ifdef DA_WITH_SYMBOLIC
        if constexpr (std::is_same_v<T, SymEngine::Expression>) {
            if (!da::is_zero(p[i])) return false;
        } else
#endif
        {
            if (std::abs(p[i]) > zero) return false;
        }
    return true;
}

// ===========================================================================
// ad_norm
// ===========================================================================
template<class T>
double ad_norm(Layout& /*layout*/, Pool<T>& pool, unsigned iv)
{
    double norm = 0.0;
    const T* p = pool.slot(iv);
    unsigned l = pool.len(iv);
    for (unsigned i = 0; i < l; ++i) {
        double v;
        if (scalar_is_near_zero(p[i])) { v = 0.0; }
        else {
#ifdef DA_WITH_SYMBOLIC
            if constexpr (std::is_same_v<T, SymEngine::Expression>) {
                v = 0.0; // norm not meaningful for purely symbolic
            } else
#endif
            { v = std::abs(scalar_to_double(p[i])); }
        }
        if (v > norm) norm = v;
    }
    return norm;
}

// ===========================================================================
// ad_weighted_norm
// ===========================================================================
template<class T>
double ad_weighted_norm(Layout& layout, Pool<T>& pool,
                        unsigned iv, double w)
{
    unsigned gnv = layout.num_vars();
    unsigned gnd = layout.max_order();
    const unsigned int* bptr = layout.base();

    std::vector<double> ww(gnd+1, 1.0);
    for (unsigned i = 1; i < gnd+1; ++i) ww[i] = ww[i-1] * w;

    const T* pv = pool.slot(iv);
    const unsigned int* p = bptr;
    unsigned l = pool.len(iv);
    double norm = 0.0;

    for (unsigned i = 0; i < l; ++i) {
        if (scalar_is_near_zero(pv[i])) {
            p += gnv;
            continue;
        }
        int order = 0;
        for (unsigned j = 0; j < gnv-1; ++j) {
            order += static_cast<int>(*p - *(p+1));
            ++p;
        }
        order += static_cast<int>(*p++);
        double value;
#ifdef DA_WITH_SYMBOLIC
        if constexpr (std::is_same_v<T, SymEngine::Expression>) { value = 0.0; }
        else
#endif
        { value = std::abs(scalar_to_double(pv[i])) * ww[static_cast<unsigned>(order)]; }
        if (value > norm) norm = value;
    }
    return norm;
}

// ===========================================================================
// Explicit instantiation for double
// ===========================================================================

#define INST(T) \
template void ad_copy<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_reset<T>(Layout&, Pool<T>&, unsigned); \
template void ad_reset_vector<T>(Layout&, Pool<T>&, unsigned); \
template void ad_const<T>(Layout&, Pool<T>&, unsigned, T); \
template void ad_add<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_sub<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_mult_const<T>(Layout&, Pool<T>&, unsigned, T); \
template void ad_add_const<T>(Layout&, Pool<T>&, unsigned, T); \
template void ad_div_c<T>(Layout&, Pool<T>&, unsigned, T); \
template void ad_mult_c<T>(Layout&, Pool<T>&, unsigned, T, unsigned); \
template void ad_mult<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned); \
template void ad_c_div<T>(Layout&, Pool<T>&, unsigned, T, unsigned); \
template void ad_div<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned); \
template void ad_pok<T>(Layout&, Pool<T>&, unsigned, const int*, std::size_t, T); \
template void ad_pek<T>(Layout&, Pool<T>&, unsigned, const int*, std::size_t, T&); \
template void ad_elem<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned*, T&); \
template void ad_var<T>(Layout&, Pool<T>&, unsigned, T, unsigned); \
template void ad_sqrt<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_exp<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_log<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_sin<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_cos<T>(Layout&, Pool<T>&, unsigned, unsigned); \
template void ad_der<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned); \
template void ad_int<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned); \
template void ad_composition<T>(Layout&, Pool<T>&, std::vector<unsigned>&, \
    std::vector<unsigned>&, std::vector<unsigned>&); \
template void ad_composition<T>(Layout&, Pool<T>&, std::vector<unsigned>&, \
    std::vector<double>&, std::vector<double>&); \
template void ad_composition<T>(Layout&, Pool<T>&, std::vector<unsigned>&, \
    std::vector<std::complex<double>>&, std::vector<std::complex<double>>&); \
template void ad_substitute<T>(Layout&, Pool<T>&, unsigned, unsigned, unsigned, unsigned); \
template void ad_substitute<T>(Layout&, Pool<T>&, std::vector<unsigned>&, \
    std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&); \
template void ad_clean<T>(Layout&, Pool<T>&, unsigned, double); \
template int ad_n_element<T>(Layout&, Pool<T>&, unsigned); \
template bool ad_zero_check<T>(Layout&, Pool<T>&, unsigned, double); \
template double ad_norm<T>(Layout&, Pool<T>&, unsigned); \
template double ad_weighted_norm<T>(Layout&, Pool<T>&, unsigned, double);

INST(double)

#ifdef DA_WITH_SYMBOLIC
INST(SymEngine::Expression)
#endif

#undef INST

} // namespace detail
} // namespace da
