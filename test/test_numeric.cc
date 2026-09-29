/**
 * @file test_numeric.cc
 * @brief Port of ref/tpsa/test/tests.cc against the new da:: API.
 *
 * Changes from the reference:
 *   - DAVector -> da::NDA (= da::DAVector<double>)
 *   - da[i]    -> da::base[i]
 *   - da_init, da_clear, etc. -> da:: namespace
 *   - Tests run from build dir, but reference data files are in the source
 *     test/ directory.  We use a relative path based on CATCH_CONFIG_SOURCE_PATH
 *     or fall back to a hard-coded path for CI.
 *   - compare_da_with_file / compare_cd_with_file use the new API.
 *   - We DO NOT define CATCH_CONFIG_MAIN here (it's in catch_main.cc).
 *
 * NOTE: The file-comparison tests require the working directory to contain
 * (or have access to) the .txt reference files.  The test executable must be
 * run from the directory containing those files, or the tests will skip.
 *
 * Adaptations:
 *   - The local read_da_from_file_c helper is replaced by da::read_da_from_file
 *     and da::compare_da_vectors.
 *   - Using namespace std::complex_literals for 1i.
 */

#include "catch.hpp"
#include "da/da.h"

#include <algorithm>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using da::NDA;
using std::complex;
using std::string;
using namespace std::complex_literals;
using std::vector;

// ---------------------------------------------------------------------------
// Helper: check that a file exists (for test skipping)
// ---------------------------------------------------------------------------
static bool file_exists(const std::string& fname)
{
    std::ifstream f(fname);
    return f.good();
}

// ---------------------------------------------------------------------------
// Local helpers for tests (matching ref tests.cc API)
// ---------------------------------------------------------------------------
static bool read_da_from_file_c(string filename, NDA& d)
{
    return da::read_da_from_file(filename, d);
}

static bool compare_da_vectors_c(NDA& a, NDA& b, double eps)
{
    return da::compare_da_vectors(a, b, eps);
}

// ---------------------------------------------------------------------------
// TEST CASES (mirror of ref/tpsa/test/tests.cc)
// ---------------------------------------------------------------------------

TEST_CASE("INITIALIZE DA ENVIRONMENT (numeric)", "[numeric_init]") {
    int da_dim   = 3;
    int da_order = 4;
    int n_vec    = 400;
    int init = da::da_init(static_cast<unsigned>(da_order),
                           static_cast<unsigned>(da_dim),
                           static_cast<unsigned>(n_vec));
    REQUIRE(init == 0);
}

TEST_CASE("DA FUNCTIONS (numeric)", "[numeric]") {
    // Requires the reference .txt files to be present in CWD.
    if (!file_exists("sqrt_da.txt")) {
        WARN("Reference data files not found in CWD -- skipping file-compare tests");
        return;
    }

    NDA x = NDA(1.0) + da::base[0] + NDA(2.0)*da::base[1] + NDA(5.0)*da::base[2];
    double eps = 1e-14;

    NDA y = da::sqrt(x);
    SECTION("SQRT") {
        REQUIRE(da::compare_da_with_file("sqrt_da.txt", y, eps));
    }

    y = da::log(x);
    SECTION("LOG") {
        REQUIRE(da::compare_da_with_file("log_da.txt", y, eps));
    }

    y = da::log(x);
    y = da::da_int(y, 1);
    NDA m;
    read_da_from_file_c("da_int.txt", m);
    SECTION("DA_INT") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }

    y = da::log(x);
    y = da::da_der(y, 1);
    read_da_from_file_c("da_der.txt", m);
    SECTION("DA_DER") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }

    y = da::pow(x, 3);
    read_da_from_file_c("pow3_da.txt", m);
    SECTION("POW(x,3)") {
        REQUIRE(compare_da_vectors_c(m, y, eps));
    }

    y = da::pow(x, 0.3);
    SECTION("POW(x,0.3)") {
        REQUIRE(da::compare_da_with_file("pow0p3_da.txt", y, eps));
    }

    y = da::exp(x);
    SECTION("EXP") {
        REQUIRE(da::compare_da_with_file("exp_da.txt", y, eps));
    }

    y = da::exp(x);
    NDA z;
    da::da_substitute(y, 0, 1.0, z);
    SECTION("SUBSTITUTE A NUMBER") {
        REQUIRE(da::compare_da_with_file("substitute_number.txt", z, eps));
    }

    y = da::exp(x);
    da::da_substitute(y, 0, x, z);
    SECTION("SUBSTITUTE A DA VECTOR") {
        REQUIRE(da::compare_da_with_file("substitute_da_vector.txt", z, eps));
    }

    y = da::exp(x);
    std::vector<NDA> lv(2);
    lv[0] = da::sin(x);
    lv[1] = da::cos(x);
    std::vector<unsigned int> idx{0, 1};
    da::da_substitute(y, idx, lv, z);
    SECTION("SUBSTITUTE MULTIPLE DA VECTORS") {
        REQUIRE(da::compare_da_with_file("substitute_multiple_da_vectors.txt", z, eps));
    }

    y = da::exp(x);
    std::vector<NDA> lx(3);
    std::vector<NDA> ly(3);
    lx[0] = x;
    lx[1] = y;
    lx[2] = da::sinh(x);

    da::da_substitute(lx, idx, lv, ly);
    SECTION("BUNCH SUBSTITUTION") {
        REQUIRE(da::compare_da_with_file("bunch_substitution_0.txt", ly[0], eps));
        REQUIRE(da::compare_da_with_file("bunch_substitution_1.txt", ly[1], eps));
        REQUIRE(da::compare_da_with_file("bunch_substitution_2.txt", ly[2], eps));
    }

    std::vector<NDA> lu(3);
    lu[0] = da::sin(x);
    lu[1] = da::cos(x);
    lu[2] = da::tan(x);
    da::da_composition(lx, lu, ly);
    SECTION("DA COMPOSITION") {
        REQUIRE(da::compare_da_with_file("da_composition_0.txt", ly[0], eps));
        REQUIRE(da::compare_da_with_file("da_composition_1.txt", ly[1], eps));
        REQUIRE(da::compare_da_with_file("da_composition_2.txt", ly[2], eps));
    }
}

TEST_CASE("CD FUNCTIONS (numeric)", "[numeric_cd]") {
    // Self-initialize so this case can also be run in isolation
    // (Catch2 runs cases independently; do not rely on an earlier case).
    da::da_init(4, 3, 400);

    if (!file_exists("cd_calculation_0.txt")) {
        WARN("Reference data files not found in CWD -- skipping CD tests");
        return;
    }

    double eps = 1e-14;
    NDA x1, x2, x3, x4;
    x1 = da::base[0] + NDA(2.0)*da::base[1] + NDA(3.0)*da::base[2];
    x2 = da::sin(x1);
    x1 = da::cos(x1);

    x3 = NDA(0.5)*da::base[0] + NDA(4.0)*da::base[1] + NDA(2.7)*da::base[2];
    x4 = da::sin(x3);
    x3 = da::cos(x3);

    auto y1 = x1 + x2 * 1i;
    auto y2 = x3 + x4 * 1i;

    SECTION("CD FUNDAMENTAL CALCULATIONS") {
        auto r = y1 + y2;
        REQUIRE(da::compare_cd_with_file("cd_calculation_0.txt", r, eps));
        r = y1 - y2;
        REQUIRE(da::compare_cd_with_file("cd_calculation_1.txt", r, eps));
        r = y1 * y2;
        REQUIRE(da::compare_cd_with_file("cd_calculation_2.txt", r, eps));
        r = y1 / y2;
        REQUIRE(da::compare_cd_with_file("cd_calculation_3.txt", r, eps));
    }

    std::vector<NDA> mmap;
    mmap.push_back(x1);
    mmap.push_back(x2);

    std::vector<complex<NDA>> cnmap;
    cnmap.push_back(y1);
    cnmap.push_back(y2);
    cnmap.push_back(y1 * y2);

    std::vector<complex<NDA>> comap(2);
    da::cd_composition(mmap, cnmap, comap);

    SECTION("DA COMPOSITION CD") {
        REQUIRE(da::compare_cd_with_file("da_composition_cd_0.txt", comap[0], eps));
        REQUIRE(da::compare_cd_with_file("da_composition_cd_1.txt", comap[1], eps));
    }

    std::vector<complex<NDA>> cmmap;
    cmmap.push_back(x1 + 1i * da::exp(x1));
    cmmap.push_back(x2 + 1i * da::exp(x2));

    da::cd_composition(cmmap, cnmap, comap);
    SECTION("CD COMPOSITION CD") {
        REQUIRE(da::compare_cd_with_file("cd_composition_cd_0.txt", comap[0], eps));
        REQUIRE(da::compare_cd_with_file("cd_composition_cd_1.txt", comap[1], eps));
    }

    mmap.push_back(x1 + NDA(0.33) * x2);
    da::cd_composition(cmmap, mmap, comap);
    SECTION("CD COMPOSITION DA") {
        REQUIRE(da::compare_cd_with_file("cd_composition_da_0.txt", comap[0], eps));
        REQUIRE(da::compare_cd_with_file("cd_composition_da_1.txt", comap[1], eps));
    }
}

// ===========================================================================
// Copying a DA vector while the order is temporarily lowered must truncate.
// ad_copy used to ignore Layout::full_len(), which made inv_map treat every
// map as purely linear.
// ===========================================================================
TEST_CASE("da_change_order truncates a copy", "[numeric_order]") {
    da::da_init(5, 2, 1000);
    // Every vector lives in an inner scope: da_clear() destroys the env, so no
    // DAVector may outlive it (see test_multienv.cc).
    {
        NDA x = da::base[0];
        NDA y = da::base[1];
        NDA v = 2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y;

        da::da_change_order(1);
        NDA linear = v;             // keeps only the linear terms
        da::da_restore_order();

        REQUIRE(NDA::order() == 5);

        NDA nonlinear = v - linear;
        REQUIRE(nonlinear.norm() > 1e-12);                  // quadratic part survives

        std::vector<int> c(2, 0);
        c[0] = 2;
        REQUIRE(std::fabs(nonlinear.element(c) - 0.1) < 1e-14);
        c[0] = 1; c[1] = 1;
        REQUIRE(std::fabs(nonlinear.element(c) - 0.05) < 1e-14);
        c[0] = 1; c[1] = 0;
        REQUIRE(std::fabs(nonlinear.element(c)) < 1e-14);   // linear part removed
    }
    da::da_clear();
}

// ===========================================================================
// inv_map: the inverse composed with the original map is the identity.
// ===========================================================================
TEST_CASE("inv_map inverts a nonlinear map", "[numeric_inv]") {
    da::da_init(5, 2, 1000);
    {
        NDA x = da::base[0];
        NDA y = da::base[1];

        std::vector<NDA> map = { 2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y,
                                -0.4 * x + 1.5 * y + 0.2 * y * y };
        std::vector<NDA> inv(2);
        da::inv_map(map, 2, inv);

        std::vector<NDA> composed(2);
        da::da_composition(map, inv, composed);
        for (int i = 0; i < 2; ++i) {
            NDA residual = composed[i] - da::base[i];
            REQUIRE(residual.norm() < 1e-12);
        }

        da::da_composition(inv, map, composed);         // other direction too
        for (int i = 0; i < 2; ++i) {
            NDA residual = composed[i] - da::base[i];
            REQUIRE(residual.norm() < 1e-12);
        }
    }
    da::da_clear();
}

// ===========================================================================
// Library functions that lower the order internally must hand back the order
// the caller had set, not the order given to da_init.
// ===========================================================================
TEST_CASE("inv_map keeps a caller-lowered order", "[numeric_order]") {
    da::da_init(5, 2, 1000);
    {
        NDA x = da::base[0];
        NDA y = da::base[1];
        std::vector<NDA> map = { 2.0 * x + 0.3 * y + 0.1 * x * x,
                                -0.4 * x + 1.5 * y + 0.2 * y * y };
        da::da_change_order(3);
        std::vector<NDA> inv(2);
        da::inv_map(map, 2, inv);
        REQUIRE(da::da_current_env().layout().max_order() == 3);
        da::da_restore_order();
    }
    da::da_clear();
}

TEST_CASE("da_composition keeps a caller-lowered order", "[numeric_order]") {
    da::da_init(5, 2, 1000);
    {
        NDA x = da::base[0];
        NDA y = da::base[1];
        std::vector<NDA> f = { x + 0.5 * x * y + 0.2 * y * y * y, y - 0.3 * x * x };
        std::vector<NDA> g = { 1.0 + x + 0.1 * y * y, 2.0 * y + 0.4 * x * y };

        std::vector<NDA> full(2);
        da::da_composition(f, g, full);       // at order 5

        da::da_change_order(3);
        std::vector<NDA> low(2);
        da::da_composition(f, g, low);
        REQUIRE(da::da_current_env().layout().max_order() == 3);

        // Expected: the order-5 result truncated to order 3.
        for (int i = 0; i < 2; ++i) {
            NDA expected = full[i];           // copy truncates at the current order
            REQUIRE(da::compare_da_vectors(expected, low[i], 1e-14));
        }
        da::da_restore_order();
    }
    da::da_clear();
}

TEST_CASE("multi-base da_substitute keeps a caller-lowered order", "[numeric_order]") {
    da::da_init(5, 2, 1000);
    {
        NDA x = da::base[0];
        NDA y = da::base[1];
        std::vector<NDA> f = { x + 0.5 * x * y + 0.2 * y * y * y, y - 0.3 * x * x };
        // Substitute base 0 only: terms that keep y take the internal
        // lower-order product path in ad_substitute.
        std::vector<unsigned int> idx{0};
        std::vector<NDA> g = { x + 0.1 * y * y + 0.3 * x * y };

        std::vector<NDA> full(2);
        da::da_substitute(f, idx, g, full);   // at order 5

        da::da_change_order(3);
        std::vector<NDA> low(2);
        da::da_substitute(f, idx, g, low);
        REQUIRE(da::da_current_env().layout().max_order() == 3);

        for (int i = 0; i < 2; ++i) {
            NDA expected = full[i];           // copy truncates at the current order
            REQUIRE(da::compare_da_vectors(expected, low[i], 1e-14));
        }
        da::da_restore_order();
    }
    da::da_clear();
}
