/**
 * @file benchmark_composition.cc
 * @brief Benchmark: single composition of 6-base DA vectors (Table 1 from ref/tpsa README).
 *
 * Reproduces Table 1: time a single composition of ONE DA vector of 6 bases
 * with 6 DA vectors, for orders 2, 4, 6 (8, 10 if pool allows).
 * Prints time for BOTH the MiraDAC engine and the reference ref/tpsa engine
 * side by side to confirm no regression.
 *
 * Build (after building /tmp/libref_tpsa.a):
 *   cmake -S . -B build -DWITH_SYMBOLIC=OFF
 *   cmake --build build --target benchmark_composition -j4
 *
 * Run:
 *   ./build/examples/benchmark_composition
 *
 * The reference lib /tmp/libref_tpsa.a must exist.  It is built from:
 *   mkdir -p /tmp/libref_build
 *   g++ -std=c++17 -O2 -c -I ref/tpsa/include \
 *       ref/tpsa/src/tpsa_extend.cc -o /tmp/libref_build/tpsa_extend.o
 *   ar rcs /tmp/libref_tpsa.a /tmp/libref_build/tpsa_extend.o
 */

// ===== Reference engine headers =====
#include "tpsa_extend.h"
// Undefine the Fortran-convention macros introduced by tpsa.h on Linux
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

// ===== MiraDAC headers =====
#include "da/da.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers for the reference engine
// ---------------------------------------------------------------------------
static void ref_init(unsigned nv, unsigned nd, unsigned pool_size)
{
    TNVND v = static_cast<TNVND>(nv);
    TNVND d = static_cast<TNVND>(nd);
    ad_init_(&v, &d);
    ad_reserve_(pool_size);
    ad_generate_order_table();
}

static unsigned comb(unsigned n, unsigned k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n-k) k = n-k;
    unsigned r = 1;
    for (unsigned i = 0; i < k; ++i) {
        r = r * (n-i) / (i+1);
    }
    return r;
}

// ---------------------------------------------------------------------------
// Benchmark one order with the MiraDAC engine
// ---------------------------------------------------------------------------
static double bench_miradac(unsigned order, unsigned nv, unsigned pool_size, unsigned reps)
{
    using namespace da;
    da_init(order, nv, pool_size);

    // Build one DA vector: identity map on first variable (uses base[0])
    std::vector<NDA> ivec(1);
    ivec[0] = base[0];  // x_0 = base[0]

    // Build 6 mapping DA vectors: each is the corresponding basis variable
    std::vector<NDA> maps(nv);
    for (unsigned i = 0; i < nv; ++i)
        maps[i] = base[i];

    std::vector<NDA> ovec(1);

    // Warm-up
    da_composition(ivec, maps, ovec);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned r = 0; r < reps; ++r) {
        da_composition(ivec, maps, ovec);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    da_clear();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    return elapsed / static_cast<double>(reps);
}

// ---------------------------------------------------------------------------
// Benchmark one order with the reference engine
// ---------------------------------------------------------------------------
static double bench_ref(unsigned order, unsigned nv, unsigned pool_size, unsigned reps)
{
    ref_init(nv, order, pool_size);

    // Allocate 1 input + nv map + 1 output vector
    std::vector<TVEC> ivec(1), maps(nv), ovec(1);
    ad_alloc_(&ivec[0]);
    for (unsigned i = 0; i < nv; ++i) ad_alloc_(&maps[i]);
    ad_alloc_(&ovec[0]);

    // ivec[0] = x_0 (first basis variable)
    {
        double x0 = 0.0;
        unsigned ibv = 1;
        ad_var_(&ivec[0], &x0, &ibv);
    }
    // maps[i] = base[i+1]
    for (unsigned i = 0; i < nv; ++i) {
        double x0 = 0.0;
        unsigned ibv = i+1;
        ad_var_(&maps[i], &x0, &ibv);
    }

    // Warm-up
    ad_composition(ivec, maps, ovec);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned r = 0; r < reps; ++r) {
        ad_composition(ivec, maps, ovec);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    // Free
    ad_free_(&ivec[0]);
    for (unsigned i = 0; i < nv; ++i) ad_free_(&maps[i]);
    ad_free_(&ovec[0]);

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    return elapsed / static_cast<double>(reps);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
    const unsigned nv = 6;
    // (order, pool_size, reps) — fewer reps for high orders
    struct Config { unsigned order; unsigned pool_size; unsigned reps; };
    std::vector<Config> configs = {
        { 2,  1000, 10000 },
        { 4,  2000,  1000 },
        { 6,  5000,   200 },
        { 8, 10000,    20 },
        {10, 20000,     5 },
    };

    std::cout << "\n=== MiraDAC vs ref/tpsa: single composition benchmark (nv=" << nv << ") ===\n";
    std::cout << std::setw(7)  << "Order"
              << std::setw(10) << "#Terms"
              << std::setw(18) << "MiraDAC (s)"
              << std::setw(18) << "ref/tpsa (s)"
              << std::setw(12) << "Ratio"
              << "\n";
    std::cout << std::string(65, '-') << "\n";

    for (auto& cfg : configs) {
        unsigned terms = comb(nv + cfg.order, cfg.order);
        double t_new = bench_miradac(cfg.order, nv, cfg.pool_size, cfg.reps);
        double t_ref = bench_ref   (cfg.order, nv, cfg.pool_size, cfg.reps);
        double ratio = t_new / t_ref;

        std::cout << std::setw(7)  << cfg.order
                  << std::setw(10) << terms
                  << std::setw(18) << std::scientific << std::setprecision(3) << t_new
                  << std::setw(18) << std::scientific << std::setprecision(3) << t_ref
                  << std::setw(12) << std::fixed      << std::setprecision(2) << ratio
                  << "\n";
    }
    std::cout << "\nRatio < 2 indicates no significant regression.\n\n";
    return 0;
}
