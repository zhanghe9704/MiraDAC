/**
 * @file test_kernels.cc
 * @brief Bit-exact comparison of da::detail kernels (double) against the reference engine.
 *
 * Build (from repo root):
 *   g++ -std=c++17 -O2 -g \
 *       -I test/ -I include/ -I ref/tpsa/include/ \
 *       test/catch_main.cc test/test_kernels.cc \
 *       /tmp/da_build/libdaStatic.a /tmp/libref_tpsa.a \
 *       -o /tmp/test_kernels
 *
 * Symbol disambiguation:
 *   tpsa.h #defines ad_copy -> ad_copy_ etc. on non-Windows.
 *   We undefine those macros right after including tpsa_extend.h so
 *   our da::detail calls are not renamed.  The underscore variants
 *   (ad_copy_, ad_mult_, etc.) are already declared by tpsa.h and
 *   present in the reference static lib.
 *
 * Reference engine lifecycle note:
 *   ad_clear() has a bug — it deletes uninitialized prdidx entries for
 *   the highest-order monomials, causing a heap corruption crash on the
 *   second call.  We therefore call ad_init_ + ad_reserve_ only ONCE per
 *   (NV, ND) configuration, and use ad_pool_clean(0) to reset the pool
 *   between sections.  ad_clear() is called only once in the final
 *   "kernels_cleanup" test case.
 */

// ===== Reference engine =====
#include "tpsa_extend.h"   // also #defines ad_copy -> ad_copy_ etc. on Linux
// Undefine the Fortran-convention macros so da::detail calls below are clean.
#undef ad_print
#undef ad_elem
#undef ad_fill_ran
#undef ad_tra
#undef ad_shift
#undef ad_save_block
#undef ad_read_block
#undef ad_nvar
#undef ad_length
#undef ad_subst
#undef ad_cos
#undef ad_sin
#undef ad_log
#undef ad_exp
#undef ad_sqrt
#undef ad_abs
#undef ad_div_c
#undef ad_c_div
#undef ad_mult_const
#undef ad_add_const
#undef ad_div
#undef ad_mult
#undef ad_sub
#undef ad_reset
#undef ad_resetvars
#undef ad_pok
#undef ad_pek
#undef ad_truncate
#undef ad_var
#undef ad_count
#undef ad_const
#undef ad_free
#undef ad_add
#undef ad_copy
#undef ad_clean
#undef ad_alloc
#undef ad_reserve
#undef ad_init

// External globals from tpsa.cpp (defined in the reference static lib)
#include <vector>
extern std::vector<double*> advec;
extern std::vector<unsigned int> adveclen;

// Convenience wrappers calling the underscore-named reference symbols.
// (tpsa.h declarations, after macro expansion, are void ad_copy_(...) etc.)
static inline void rref_alloc(unsigned& i)         { ad_alloc_(&i); }
static inline void rref_free(unsigned i)            { ad_free_(&i); }
static inline void rref_copy(unsigned s, unsigned d){ ad_copy_(&s, &d); }
static inline void rref_reset(unsigned iv)          { ad_reset_(&iv); }
static inline void rref_const(unsigned iv, double r){ ad_const_(&iv, &r); }
static inline void rref_add(unsigned d, unsigned s) { ad_add_(&d, &s); }
static inline void rref_sub(unsigned d, unsigned s) { ad_sub_(&d, &s); }
static inline void rref_mult(unsigned l, unsigned r, unsigned d) { ad_mult_(&l, &r, &d); }
static inline void rref_mult_const(unsigned iv, double c) { ad_mult_const_(&iv, &c); }
static inline void rref_div_c(unsigned iv, double c)      { ad_div_c_(&iv, &c); }
static inline void rref_c_div(unsigned iv, double c, unsigned ivret) { ad_c_div_(&iv, &c, &ivret); }
static inline void rref_div(unsigned l, unsigned r, unsigned d) { ad_div_(&l, &r, &d); }
static inline void rref_pok(unsigned iv, int* c, size_t n, double x) { ad_pok_(&iv, c, &n, &x); }
static inline void rref_pek(unsigned iv, int* c, size_t n, double& x) { ad_pek_(&iv, c, &n, &x); }
static inline void rref_var(unsigned iv, double x0, unsigned ibv) { ad_var_(&iv, &x0, &ibv); }
static inline void rref_sqrt(unsigned iv, unsigned iret) { ad_sqrt_(&iv, &iret); }
static inline void rref_exp(unsigned iv, unsigned iret)  { ad_exp_(&iv, &iret); }
static inline void rref_log(unsigned iv, unsigned iret)  { ad_log_(&iv, &iret); }
static inline void rref_sin(unsigned iv, unsigned iret)  { ad_sin_(&iv, &iret); }
static inline void rref_cos(unsigned iv, unsigned iret)  { ad_cos_(&iv, &iret); }

// ===== New engine =====
#include "da/layout.h"
#include "da/pool.h"
#include "da/env.h"
#include "da/engine.h"

// ===== Catch2 =====
#include "catch.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>
#include <random>
#include <sstream>
#include <iomanip>

static constexpr double EXACT_TOL = 0.0;
static constexpr double TRANS_TOL = 1e-12;

// Global state — persists across sections within one TEST_CASE run.
static da::DAEnv* g_env = nullptr;
static unsigned g_nv_init = 0;
static unsigned g_nd_init = 0;
static unsigned g_pool_size_init = 0;

/**
 * ref_setup: initialise the reference engine ONLY if (nv,nd) differ from
 * the current configuration.  Otherwise just flush the pool with
 * ad_pool_clean(0) so that section re-runs start from a clean state without
 * calling ad_init_ again (which would trigger the prdidx bug in ad_clear).
 */
static void ref_setup(unsigned nv, unsigned nd, unsigned pool_size = 2000)
{
    // The reference engine has two bugs we must work around:
    //
    //  1. ad_clear() deletes prdidx rows that were never initialised
    //     (for monomials of the highest order), causing free(): invalid size
    //     on the second call.  We therefore NEVER call ad_clear() during
    //     tests — only in kernels_cleanup after all tests finish.
    //
    //  2. ad_pool_clean(0) has an unsigned underflow bug (size_t i = idx-1
    //     wraps to SIZE_MAX when idx==0), so it never rebuilds adlist.
    //     We work around it by calling ad_reserve_ instead.
    //
    // ad_init_ is safe to call repeatedly: it writes advecpool = NULL
    // (preventing double-free in ad_reserve_), then allocates base/prdidx
    // fresh (the old pointers are leaked — acceptable in tests).

    if (g_nv_init != nv || g_nd_init != nd || g_pool_size_init != pool_size) {
        // Configuration change or first call: reinitialise everything.
        TNVND v = static_cast<TNVND>(nv);
        TNVND d = static_cast<TNVND>(nd);
        ad_init_(&v, &d);
        ad_reserve_(pool_size);
        ad_generate_order_table();
        g_nv_init        = nv;
        g_nd_init        = nd;
        g_pool_size_init = pool_size;
    } else {
        // Same config: re-run ad_reserve_ to rebuild the free-list cleanly.
        ad_reserve_(pool_size);
    }
}

/**
 * new_setup: (re)create the da::DAEnv.  Called every section run to give a
 * fresh Pool<double> with all slots free.
 */
static void new_setup(unsigned nv, unsigned nd, unsigned pool_size = 2000)
{
    delete g_env;
    g_env = new da::DAEnv(nd, nv, pool_size, true);
}

static void fill_random(unsigned rslot, unsigned nslot, da::DAEnv* env,
                         std::mt19937& rng, double scale = 1.0,
                         bool first_positive = true)
{
    std::uniform_real_distribution<double> dist(-scale, scale);
    unsigned fl = env->layout().full_len();
    double* rp = advec[rslot];
    double* np = env->pool<double>().slot(nslot);
    for (unsigned i = 0; i < fl; ++i) { double v = dist(rng); rp[i] = v; np[i] = v; }
    if (first_positive && fl > 0) {
        double v = std::abs(dist(rng)) + 0.5;
        rp[0] = v; np[0] = v;
    }
    adveclen[rslot] = fl;
    env->pool<double>().set_len(nslot, fl);
}

static double compare_slots(unsigned rslot, unsigned nslot, da::DAEnv* env,
                              double tol, const std::string& label)
{
    unsigned fl = env->layout().full_len();
    const double* rp = advec[rslot];
    const double* np = env->pool<double>().slot(nslot);
    double maxdiff = 0.0;
    for (unsigned i = 0; i < fl; ++i) {
        double d = std::abs(rp[i] - np[i]);
        if (d > maxdiff) maxdiff = d;
    }
    INFO(label << " max_diff=" << std::scientific << maxdiff);
    if (tol == EXACT_TOL) REQUIRE(maxdiff == 0.0);
    else                   REQUIRE(maxdiff <= tol);
    return maxdiff;
}

static std::pair<unsigned,unsigned> both_alloc(da::DAEnv* env)
{
    unsigned rs, ns;
    rref_alloc(rs);
    ns = env->pool<double>().alloc();
    return {rs, ns};
}
static void both_free(da::DAEnv* env, unsigned rs, unsigned ns)
{
    rref_free(rs);
    env->pool<double>().free(ns);
}

// ===========================================================================
// (3,4): nv=3, nd=4, full_len=35
// ===========================================================================
TEST_CASE("kernels_3_4: nv=3 nd=4", "[kernels]")
{
    unsigned NV = 3, ND = 4, POOL = 2000;
    ref_setup(NV, ND, POOL);
    new_setup(NV, ND, POOL);
    da::Layout& L = g_env->layout();
    da::Pool<double>& P = g_env->pool<double>();
    REQUIRE(L.full_len() == 35u);
    std::mt19937 rng(42);

    SECTION("copy") {
        auto [rs, ns] = both_alloc(g_env);
        auto [rd, nd2] = both_alloc(g_env);
        fill_random(rs, ns, g_env, rng);
        rref_copy(rs, rd);
        da::detail::ad_copy(L, P, ns, nd2);
        compare_slots(rd, nd2, g_env, EXACT_TOL, "ad_copy");
        both_free(g_env, rs, ns); both_free(g_env, rd, nd2);
    }
    SECTION("reset") {
        auto [rs, ns] = both_alloc(g_env);
        fill_random(rs, ns, g_env, rng);
        rref_reset(rs); da::detail::ad_reset(L, P, ns);
        compare_slots(rs, ns, g_env, EXACT_TOL, "ad_reset");
        both_free(g_env, rs, ns);
    }
    SECTION("const") {
        auto [rs, ns] = both_alloc(g_env);
        fill_random(rs, ns, g_env, rng);
        double c = 3.14;
        rref_const(rs, c); da::detail::ad_const(L, P, ns, c);
        compare_slots(rs, ns, g_env, EXACT_TOL, "ad_const");
        both_free(g_env, rs, ns);
    }
    SECTION("add") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        fill_random(rb, nb, g_env, rng, 1.0, false);
        rref_add(ra, rb); da::detail::ad_add(L, P, na, nb);
        compare_slots(ra, na, g_env, EXACT_TOL, "ad_add");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("sub") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        fill_random(rb, nb, g_env, rng, 1.0, false);
        rref_sub(ra, rb); da::detail::ad_sub(L, P, na, nb);
        compare_slots(ra, na, g_env, EXACT_TOL, "ad_sub");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("mult_const") {
        auto [ra, na] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        double c = 2.718;
        rref_mult_const(ra, c); da::detail::ad_mult_const(L, P, na, c);
        compare_slots(ra, na, g_env, EXACT_TOL, "ad_mult_const");
        both_free(g_env, ra, na);
    }
    SECTION("mult") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        auto [rc, nc] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        fill_random(rb, nb, g_env, rng, 1.0, false);
        rref_mult(ra, rb, rc); da::detail::ad_mult(L, P, na, nb, nc);
        compare_slots(rc, nc, g_env, EXACT_TOL, "ad_mult");
        both_free(g_env, ra, na); both_free(g_env, rb, nb); both_free(g_env, rc, nc);
    }
    SECTION("div_c") {
        auto [ra, na] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        double c = 1.5;
        rref_div_c(ra, c); da::detail::ad_div_c(L, P, na, c);
        compare_slots(ra, na, g_env, EXACT_TOL, "ad_div_c");
        both_free(g_env, ra, na);
    }
    SECTION("c_div") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, true);
        double c = 1.0;
        rref_c_div(ra, c, rb); da::detail::ad_c_div(L, P, na, c, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_c_div");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("div") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        auto [rc, nc] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, true);
        fill_random(rb, nb, g_env, rng, 0.3, true);
        rref_div(ra, rb, rc); da::detail::ad_div(L, P, na, nb, nc);
        compare_slots(rc, nc, g_env, TRANS_TOL, "ad_div");
        both_free(g_env, ra, na); both_free(g_env, rb, nb); both_free(g_env, rc, nc);
    }
    SECTION("var") {
        for (unsigned ibv = 0; ibv < NV; ++ibv) {
            auto [rs, ns] = both_alloc(g_env);
            double x0 = 1.0 + ibv * 0.5;
            rref_var(rs, x0, ibv); da::detail::ad_var(L, P, ns, x0, ibv);
            std::ostringstream ss; ss << "ad_var(ibv=" << ibv << ")";
            compare_slots(rs, ns, g_env, EXACT_TOL, ss.str());
            both_free(g_env, rs, ns);
        }
    }
    SECTION("pok/pek") {
        auto [ra, na] = both_alloc(g_env);
        int exponents[3] = {2, 1, 0};
        size_t ne = 3; double val = 7.77;
        rref_pok(ra, exponents, ne, val);
        da::detail::ad_pok(L, P, na, exponents, ne, val);
        double rval = 0, nval = 0;
        rref_pek(ra, exponents, ne, rval);
        da::detail::ad_pek(L, P, na, exponents, ne, nval);
        REQUIRE(rval == Approx(nval).epsilon(0));
        REQUIRE(nval == Approx(val).epsilon(0));
        both_free(g_env, ra, na);
    }
    SECTION("sqrt") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, false);
        advec[ra][0] = 2.0; P.slot(na)[0] = 2.0;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_sqrt(ra, rb); da::detail::ad_sqrt(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_sqrt");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("exp") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, false);
        advec[ra][0] = 0.5; P.slot(na)[0] = 0.5;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_exp(ra, rb); da::detail::ad_exp(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_exp");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("log") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, false);
        advec[ra][0] = 2.0; P.slot(na)[0] = 2.0;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_log(ra, rb); da::detail::ad_log(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_log");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("sin") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, false);
        advec[ra][0] = 0.7; P.slot(na)[0] = 0.7;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_sin(ra, rb); da::detail::ad_sin(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_sin");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("cos") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.3, false);
        advec[ra][0] = 0.7; P.slot(na)[0] = 0.7;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_cos(ra, rb); da::detail::ad_cos(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_cos");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("der") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        unsigned expo = 1;
        ad_der(&ra, &expo, &rb);
        da::detail::ad_der(L, P, na, expo, nb);
        compare_slots(rb, nb, g_env, EXACT_TOL, "ad_der");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("int") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        unsigned base_id = 0;
        ad_int(ra, base_id, rb);
        da::detail::ad_int(L, P, na, base_id, nb);
        compare_slots(rb, nb, g_env, EXACT_TOL, "ad_int");
        both_free(g_env, ra, na); both_free(g_env, rb, nb);
    }
    SECTION("norm") {
        auto [ra, na] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        double rnorm = ad_norm(ra);
        double nnorm = da::detail::ad_norm(L, P, na);
        REQUIRE(rnorm == Approx(nnorm).epsilon(1e-15));
        both_free(g_env, ra, na);
    }
    SECTION("n_element") {
        auto [ra, na] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        int rn = ad_n_element(ra);
        int nn = da::detail::ad_n_element(L, P, na);
        REQUIRE(rn == nn);
        both_free(g_env, ra, na);
    }
    SECTION("composition TPS->TPS") {
        std::vector<unsigned> riv(NV), niv(NV);
        for (unsigned i = 0; i < NV; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            riv[i] = rr; niv[i] = nn;
            fill_random(rr, nn, g_env, rng, 0.1, false);
            advec[rr][0] = 0.0; P.slot(nn)[0] = 0.0;
            advec[rr][i+1] = 1.0; P.slot(nn)[i+1] = 1.0;
            adveclen[rr] = L.full_len(); P.set_len(nn, L.full_len());
        }
        std::vector<unsigned> riv2(2), niv2(2);
        for (int i = 0; i < 2; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            riv2[i] = rr; niv2[i] = nn;
            fill_random(rr, nn, g_env, rng, 0.2, false);
            adveclen[rr] = L.full_len(); P.set_len(nn, L.full_len());
        }
        std::vector<unsigned> rov(2), nov(2);
        for (int i = 0; i < 2; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            rov[i] = rr; nov[i] = nn;
        }
        ad_composition(riv2, riv, rov);
        da::detail::ad_composition(L, P, niv2, niv, nov);
        for (int i = 0; i < 2; ++i) {
            std::ostringstream ss; ss << "composition TPS[" << i << "]";
            compare_slots(rov[i], nov[i], g_env, TRANS_TOL, ss.str());
        }
        for (unsigned i = 0; i < NV; ++i) both_free(g_env, riv[i], niv[i]);
        for (int i = 0; i < 2; ++i) {
            both_free(g_env, riv2[i], niv2[i]);
            both_free(g_env, rov[i], nov[i]);
        }
    }
    SECTION("composition TPS->double") {
        std::vector<unsigned> riv(2), niv(2);
        for (int i = 0; i < 2; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            riv[i] = rr; niv[i] = nn;
            fill_random(rr, nn, g_env, rng, 0.3, false);
            adveclen[rr] = L.full_len(); P.set_len(nn, L.full_len());
        }
        std::vector<double> dv = {0.5, -0.3, 0.1};
        std::vector<double> rov(2, 0.0), nov(2, 0.0);
        ad_composition(riv, dv, rov);
        da::detail::ad_composition(L, P, niv, dv, nov);
        for (int i = 0; i < 2; ++i)
            REQUIRE(std::abs(rov[i] - nov[i]) <= TRANS_TOL);
        for (int i = 0; i < 2; ++i) both_free(g_env, riv[i], niv[i]);
    }
    SECTION("substitute single-base") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        auto [rc, nc] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.2, false);
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        fill_random(rb, nb, g_env, rng, 0.2, false);
        advec[rb][0] = 0.0; P.slot(nb)[0] = 0.0;
        advec[rb][1] = 1.0; P.slot(nb)[1] = 1.0;
        adveclen[rb] = L.full_len(); P.set_len(nb, L.full_len());
        unsigned base_id = 0;
        ad_substitute(ra, base_id, rb, rc);
        da::detail::ad_substitute(L, P, na, base_id, nb, nc);
        compare_slots(rc, nc, g_env, TRANS_TOL, "ad_substitute single");
        both_free(g_env, ra, na); both_free(g_env, rb, nb); both_free(g_env, rc, nc);
    }
}

// ===========================================================================
// (6,6): nv=6, nd=6, full_len=924
// ===========================================================================
TEST_CASE("kernels_6_6: nv=6 nd=6 full_len=924", "[kernels]")
{
    unsigned NV = 6, ND = 6, POOL = 4000;
    ref_setup(NV, ND, POOL);
    new_setup(NV, ND, POOL);
    da::Layout& L = g_env->layout();
    da::Pool<double>& P = g_env->pool<double>();
    REQUIRE(L.full_len() == 924u);
    std::mt19937 rng(123);

    SECTION("mult (6,6)") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        auto [rc, nc] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        fill_random(rb, nb, g_env, rng, 1.0, false);
        rref_mult(ra, rb, rc); da::detail::ad_mult(L, P, na, nb, nc);
        compare_slots(rc, nc, g_env, EXACT_TOL, "ad_mult (6,6)");
        both_free(g_env, ra, na); both_free(g_env, rb, nb); both_free(g_env, rc, nc);
    }
    SECTION("composition TPS->TPS (6,6)") {
        std::vector<unsigned> riv(NV), niv(NV);
        for (unsigned i = 0; i < NV; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            riv[i] = rr; niv[i] = nn;
            fill_random(rr, nn, g_env, rng, 0.05, false);
            advec[rr][0] = 0.0; P.slot(nn)[0] = 0.0;
            advec[rr][i+1] = 1.0; P.slot(nn)[i+1] = 1.0;
            adveclen[rr] = L.full_len(); P.set_len(nn, L.full_len());
        }
        std::vector<unsigned> riv2(2), niv2(2);
        for (int i = 0; i < 2; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            riv2[i] = rr; niv2[i] = nn;
            fill_random(rr, nn, g_env, rng, 0.1, false);
            adveclen[rr] = L.full_len(); P.set_len(nn, L.full_len());
        }
        std::vector<unsigned> rov(2), nov(2);
        for (int i = 0; i < 2; ++i) {
            auto [rr, nn] = both_alloc(g_env);
            rov[i] = rr; nov[i] = nn;
        }
        ad_composition(riv2, riv, rov);
        da::detail::ad_composition(L, P, niv2, niv, nov);
        for (int i = 0; i < 2; ++i) {
            std::ostringstream ss; ss << "composition (6,6)[" << i << "]";
            compare_slots(rov[i], nov[i], g_env, TRANS_TOL, ss.str());
        }
        for (unsigned i = 0; i < NV; ++i) both_free(g_env, riv[i], niv[i]);
        for (int i = 0; i < 2; ++i) {
            both_free(g_env, riv2[i], niv2[i]);
            both_free(g_env, rov[i], nov[i]);
        }
    }
    SECTION("sqrt/exp/sin/cos (6,6)") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng, 0.1, false);
        advec[ra][0] = 1.5; P.slot(na)[0] = 1.5;
        adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
        rref_sqrt(ra, rb); da::detail::ad_sqrt(L, P, na, nb);
        compare_slots(rb, nb, g_env, TRANS_TOL, "ad_sqrt (6,6)");
        both_free(g_env, rb, nb);

        {
            auto [rb2, nb2] = both_alloc(g_env);
            fill_random(ra, na, g_env, rng, 0.1, false);
            advec[ra][0] = 0.3; P.slot(na)[0] = 0.3;
            adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
            rref_exp(ra, rb2); da::detail::ad_exp(L, P, na, nb2);
            compare_slots(rb2, nb2, g_env, TRANS_TOL, "ad_exp (6,6)");
            both_free(g_env, rb2, nb2);
        }
        {
            auto [rb2, nb2] = both_alloc(g_env);
            fill_random(ra, na, g_env, rng, 0.1, false);
            advec[ra][0] = 0.5; P.slot(na)[0] = 0.5;
            adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
            rref_sin(ra, rb2); da::detail::ad_sin(L, P, na, nb2);
            compare_slots(rb2, nb2, g_env, TRANS_TOL, "ad_sin (6,6)");
            both_free(g_env, rb2, nb2);
        }
        {
            auto [rb2, nb2] = both_alloc(g_env);
            fill_random(ra, na, g_env, rng, 0.1, false);
            advec[ra][0] = 0.5; P.slot(na)[0] = 0.5;
            adveclen[ra] = L.full_len(); P.set_len(na, L.full_len());
            rref_cos(ra, rb2); da::detail::ad_cos(L, P, na, nb2);
            compare_slots(rb2, nb2, g_env, TRANS_TOL, "ad_cos (6,6)");
            both_free(g_env, rb2, nb2);
        }
        both_free(g_env, ra, na);
    }
    SECTION("der/int round-trip (6,6)") {
        auto [ra, na] = both_alloc(g_env);
        auto [rb, nb] = both_alloc(g_env);
        auto [rc, nc] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        unsigned base_id = 2;
        ad_int(ra, base_id, rb);
        ad_der(&rb, &base_id, &rc);
        da::detail::ad_int(L, P, na, base_id, nb);
        da::detail::ad_der(L, P, nb, base_id, nc);
        compare_slots(rc, nc, g_env, EXACT_TOL, "ad_int->der round-trip");
        both_free(g_env, ra, na); both_free(g_env, rb, nb); both_free(g_env, rc, nc);
    }
    SECTION("norm/n_element/zero_check (6,6)") {
        auto [ra, na] = both_alloc(g_env);
        fill_random(ra, na, g_env, rng);
        double rnorm = ad_norm(ra);
        double nnorm = da::detail::ad_norm(L, P, na);
        REQUIRE(std::abs(rnorm - nnorm) <= 1e-15 * std::max(rnorm, 1.0));
        REQUIRE(ad_n_element(ra) == da::detail::ad_n_element(L, P, na));
        REQUIRE(ad_zero_check(ra, 1e-300) == da::detail::ad_zero_check(L, P, na, 1e-300));
        both_free(g_env, ra, na);
    }
}

// ===========================================================================
// Final cleanup.
// ad_clear() has a known bug (deletes uninitialized prdidx[order_index[nd]..
// FULL_VEC_LEN-1] entries), so we skip it and let the OS reclaim memory.
// ===========================================================================
TEST_CASE("kernels_cleanup", "[kernels]")
{
    // Intentionally do NOT call ad_clear() — it crashes due to uninitialized
    // prdidx rows in the last-order tier (ref-engine bug, read-only source).
    g_nv_init        = 0;
    g_nd_init        = 0;
    g_pool_size_init = 0;
    delete g_env;
    g_env = nullptr;
    SUCCEED("cleanup done (ad_clear skipped — ref-engine prdidx bug)");
}
