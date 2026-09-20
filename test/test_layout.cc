/**
 * @file test_layout.cc
 * @brief Unit tests for Layout (Stage 1).
 */

#include "catch.hpp"
#include "da/layout.h"
#include "da/monomial_scheme.h"

#include <vector>
#include <numeric>

using namespace da;

// -----------------------------------------------------------------------
// Helper: build a Layout and generate the order table.
// -----------------------------------------------------------------------
static Layout make_layout(unsigned int nv, unsigned int nd)
{
    UniformOrder scheme(nv, nd);
    Layout lay(scheme, /*build_order_table=*/true);
    return lay;
}

// -----------------------------------------------------------------------
// full_len checks
// -----------------------------------------------------------------------
TEST_CASE("UniformOrder full_len", "[layout]")
{
    REQUIRE(UniformOrder(3, 4).full_len() == 35);
    REQUIRE(UniformOrder(6, 6).full_len() == 924);
    REQUIRE(UniformOrder(6, 10).full_len() == 8008);

    // extra sanity: C(nv+nd, nd) for small cases
    REQUIRE(UniformOrder(1, 1).full_len() == 2);   // {1, x}
    REQUIRE(UniformOrder(2, 2).full_len() == 6);   // {1, x, y, x2, xy, y2}
    REQUIRE(UniformOrder(3, 3).full_len() == 20);
}

// -----------------------------------------------------------------------
// order_index is monotonically non-decreasing and order_index[nd+1] == full_len
// -----------------------------------------------------------------------
TEST_CASE("order_index properties", "[layout]")
{
    auto check_order_index = [](unsigned int nv, unsigned int nd) {
        UniformOrder scheme(nv, nd);
        Layout lay(scheme, false);
        const unsigned int* oi = lay.order_index();
        unsigned int fvl = lay.full_len();

        // non-decreasing
        for (unsigned int i = 0; i <= nd; ++i) {
            INFO("nv=" << nv << " nd=" << nd << " i=" << i);
            REQUIRE(oi[i] <= oi[i + 1]);
        }

        // last entry equals full_len
        REQUIRE(oi[nd + 1] == fvl);

        // strictly increasing (each order adds at least 1 element when nv>0)
        if (nv > 0 && nd > 0) {
            for (unsigned int i = 0; i <= nd; ++i) {
                INFO("nv=" << nv << " nd=" << nd << " i=" << i);
                REQUIRE(oi[i] < oi[i + 1]);
            }
        }
    };

    check_order_index(3, 4);
    check_order_index(6, 6);
    check_order_index(2, 2);
    check_order_index(1, 5);
}

// -----------------------------------------------------------------------
// Round-trip: find_index(orders(i)) == i for every i in [0, full_len)
// -----------------------------------------------------------------------
TEST_CASE("orders / find_index round-trip", "[layout]")
{
    auto check_roundtrip = [](unsigned int nv, unsigned int nd) {
        Layout lay = make_layout(nv, nd);
        unsigned int fvl = lay.full_len();

        for (unsigned int i = 0; i < fvl; ++i) {
            auto ord = lay.orders(i);
            int idx = lay.find_index(ord);
            INFO("nv=" << nv << " nd=" << nd << " i=" << i);
            REQUIRE(idx == (int)i);
        }
    };

    check_roundtrip(2, 2);
    check_roundtrip(3, 3);
    check_roundtrip(3, 4);
    check_roundtrip(6, 4);  // moderate size
}

// -----------------------------------------------------------------------
// prdidx spot-check
//
// For nv=2, nd=4:
//   Index 0: (0,0) = 1       (constant)
//   Index 1: (1,0) = x
//   Index 2: (0,1) = y
//   Index 3: (2,0) = x^2
//   Index 4: (1,1) = x*y
//   Index 5: (0,2) = y^2
//
// prdidx[1][2] should be the index of x * y = x^1*y^1 = index 4.
// prdidx[1][1] should be the index of x * x = x^2      = index 3.
// prdidx[2][2] should be the index of y * y = y^2      = index 5.
// -----------------------------------------------------------------------
TEST_CASE("prdidx spot-check", "[layout]")
{
    Layout lay = make_layout(2, 4);

    // Verify index assignments using find_index
    std::vector<int> e_const = {0, 0};
    std::vector<int> e_x     = {1, 0};
    std::vector<int> e_y     = {0, 1};
    std::vector<int> e_x2    = {2, 0};
    std::vector<int> e_xy    = {1, 1};
    std::vector<int> e_y2    = {0, 2};

    int idx_const = lay.find_index(e_const);  // 0
    int idx_x     = lay.find_index(e_x);      // 1
    int idx_y     = lay.find_index(e_y);      // 2
    int idx_x2    = lay.find_index(e_x2);     // 3
    int idx_xy    = lay.find_index(e_xy);     // 4
    int idx_y2    = lay.find_index(e_y2);     // 5

    REQUIRE(idx_const == 0);
    REQUIRE(idx_x  == 1);
    REQUIRE(idx_y  == 2);
    REQUIRE(idx_x2 == 3);
    REQUIRE(idx_xy == 4);
    REQUIRE(idx_y2 == 5);

    unsigned int** prd = lay.prdidx();

    // x * y = xy
    REQUIRE(prd[idx_x][idx_y] == (unsigned int)idx_xy);

    // x * x = x^2
    REQUIRE(prd[idx_x][idx_x] == (unsigned int)idx_x2);

    // y * y = y^2
    REQUIRE(prd[idx_y][idx_y] == (unsigned int)idx_y2);

    // Verify using find_index on summed exponent vectors
    // Product of monomial i and j: sum their exponent vectors, find the index.
    auto check_prd = [&](int i, int j) {
        auto oi = lay.orders((unsigned int)i);
        auto oj = lay.orders((unsigned int)j);
        std::vector<int> sum(oi.size());
        for (size_t k = 0; k < oi.size(); ++k) sum[k] = oi[k] + oj[k];
        int expected = lay.find_index(sum);
        INFO("prdidx[" << i << "][" << j << "] expected=" << expected
             << " got=" << prd[i][j]);
        REQUIRE(prd[i][j] == (unsigned int)expected);
    };

    // Run over all valid (i,j) pairs for nd=4, nv=2
    unsigned int nd = lay.max_order();
    const unsigned int* oi = lay.order_index();
    for (unsigned int i = 1; i < oi[nd]; ++i) {
        auto ordi = lay.orders(i);
        unsigned int ordi_total = (unsigned int)ordi[0] + (unsigned int)ordi[1];
        unsigned int M = oi[nd - ordi_total + 1];
        for (unsigned int j = 1; j < M; ++j) {
            check_prd((int)i, (int)j);
        }
    }
}

// -----------------------------------------------------------------------
// change_order / restore_order
// -----------------------------------------------------------------------
TEST_CASE("change_order and restore_order", "[layout]")
{
    UniformOrder scheme(3, 4);
    Layout lay(scheme, false);

    REQUIRE(lay.max_order() == 4);
    REQUIRE(lay.full_len()  == 35);

    int r = lay.change_order(2);
    REQUIRE(r == 0);
    REQUIRE(lay.max_order()  == 2);
    REQUIRE(lay.full_len()   == UniformOrder(3, 2).full_len());  // C(5,2)=10

    lay.restore_order();
    REQUIRE(lay.max_order() == 4);
    REQUIRE(lay.full_len()  == 35);
}
