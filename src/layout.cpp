/**
 * @file layout.cpp
 * @brief Implementation of Layout and UniformOrder.
 *
 * Algorithms ported from:
 *   ref/tpsa/src/tpsa.cpp         — gcd, comb_num, init_order_index,
 *                                    init_prod_index, init_base, choose,
 *                                    within_limit
 *   ref/tpsa/src/tpsa_extend.cc   — ADOrderTable::generate_order_table,
 *                                    ad_change_order, ad_restore_order
 *
 * The file-static globals (gnv, gnd, FULL_VEC_LEN, order_index, base,
 * prdidx, H) are now Layout members.  All algorithms are kept identical
 * to the reference; only variable access is re-wired to members.
 */

#include "da/layout.h"
#include "da/monomial_scheme.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace da {

// ======================================================================== //
// UniformOrder
// ======================================================================== //

/// Ported from comb_num() in tpsa.cpp
static unsigned int comb_num_impl(unsigned int n, unsigned int r)
{
    const unsigned int MAX_N_BASE = std::numeric_limits<unsigned int>::max() / 2;
    const unsigned int MAX_N_R = 100;

    if (n == 0 || r == 0) return 1;
    if (n < r) return 0;
    if (r > MAX_N_R && (n - r) > MAX_N_R) return 0;

    unsigned int k = (r > n - r ? n - r : r);

    std::vector<int> numerator(k), denominator(k);
    for (unsigned int i = 0; i < k; ++i) {
        numerator[i]   = (int)(n - i);
        denominator[i] = (int)(k - i);
    }

    // GCD-reduction
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            unsigned int a = (unsigned int)denominator[i];
            unsigned int b = (unsigned int)numerator[j];
            // gcd
            unsigned int aa = a, bb = b, t;
            while (bb != 0u) { t = bb; bb = aa % bb; aa = t; }
            int c = (int)aa;
            if (c > 1) {
                numerator[j]   /= c;
                denominator[i] /= c;
            }
        }
    }

    // check denominator is all-1
    for (size_t i = 0; i < k; ++i)
        if (denominator[i] != 1) return 0;

    size_t N = 1;
    for (size_t i = 0; i < k; ++i) {
        if ((MAX_N_BASE * 1.0) / (double)numerator[i] < (double)N) return 0;
        N *= (size_t)numerator[i];
    }
    return (unsigned int)N;
}

UniformOrder::UniformOrder(unsigned int nv, unsigned int nd)
    : nv_(nv), nd_(nd)
{}

unsigned int UniformOrder::full_len() const
{
    return comb_num_impl(nv_ + nd_, nd_);
}

// ======================================================================== //
// Layout — helpers
// ======================================================================== //

/*static*/ unsigned int Layout::gcd_impl(unsigned int a, unsigned int b)
{
    unsigned int t;
    while (b != 0u) { t = b; b = a % b; a = t; }
    return a;
}

/*static*/ unsigned int Layout::comb_num(unsigned int n, unsigned int r)
{
    return comb_num_impl(n, r);
}

// ======================================================================== //
// Layout — construction / destruction
// ======================================================================== //

Layout::Layout(const MonomialScheme& scheme, bool build_order_table)
    : gnv_(scheme.num_vars())
    , gnd_(scheme.max_total_order())
    , gnd_original_(scheme.max_total_order())
    , gnd_record_(0)
    , FULL_VEC_LEN_(scheme.full_len())
    , order_index_(nullptr)
    , base_(nullptr)
    , prdidx_(nullptr)
    , prdidx_rows_(0)
    , tblsize_(0)
    , order_table_valid_(false)
{
    init_order_index();
    init_base();
    init_prod_index();

    if (build_order_table)
        generate_order_table();
}

Layout::~Layout()
{
    free_tables();
}

void Layout::free_tables()
{
    if (prdidx_) {
        // Free against the allocated row count, not FULL_VEC_LEN_, which
        // change_order() may have lowered since construction.
        for (unsigned int i = 0; i < prdidx_rows_; ++i)
            delete[] prdidx_[i];
        delete[] prdidx_;
        prdidx_ = nullptr;
        prdidx_rows_ = 0;
    }
    delete[] base_;
    base_ = nullptr;
    delete[] order_index_;
    order_index_ = nullptr;
}

Layout::Layout(Layout&& o) noexcept
    : gnv_(o.gnv_)
    , gnd_(o.gnd_)
    , gnd_original_(o.gnd_original_)
    , gnd_record_(o.gnd_record_)
    , FULL_VEC_LEN_(o.FULL_VEC_LEN_)
    , order_index_(o.order_index_)
    , base_(o.base_)
    , prdidx_(o.prdidx_)
    , prdidx_rows_(o.prdidx_rows_)
    , tblsize_(o.tblsize_)
    , H_(std::move(o.H_))
    , order_table_valid_(o.order_table_valid_)
    , order_table_(std::move(o.order_table_))
    , order_index_map_(std::move(o.order_index_map_))
{
    o.order_index_ = nullptr;
    o.base_        = nullptr;
    o.prdidx_      = nullptr;
    o.prdidx_rows_ = 0;
    o.FULL_VEC_LEN_ = 0;
}

Layout& Layout::operator=(Layout&& o) noexcept
{
    if (this != &o) {
        free_tables();
        gnv_               = o.gnv_;
        gnd_               = o.gnd_;
        gnd_original_      = o.gnd_original_;
        gnd_record_        = o.gnd_record_;
        FULL_VEC_LEN_      = o.FULL_VEC_LEN_;
        order_index_       = o.order_index_;
        base_              = o.base_;
        prdidx_            = o.prdidx_;
        prdidx_rows_       = o.prdidx_rows_;
        tblsize_           = o.tblsize_;
        H_                 = std::move(o.H_);
        order_table_valid_ = o.order_table_valid_;
        order_table_       = std::move(o.order_table_);
        order_index_map_   = std::move(o.order_index_map_);

        o.order_index_ = nullptr;
        o.base_        = nullptr;
        o.prdidx_      = nullptr;
        o.prdidx_rows_ = 0;
        o.FULL_VEC_LEN_ = 0;
    }
    return *this;
}

// ======================================================================== //
// Layout — init_order_index
//   Ported from init_order_index() in tpsa.cpp (algorithm unchanged).
// ======================================================================== //

void Layout::init_order_index()
{
    unsigned int nv = gnv_;
    unsigned int nd = gnd_;

    order_index_ = new unsigned int[nd + 2];
    size_t i = 0, N = 0;
    order_index_[0] = 0;

    while (i < nd + 1) {
        N += comb_num(nv + (unsigned int)i - 1, (unsigned int)i);
        order_index_[++i] = (unsigned int)N;
    }
    order_index_[nd + 1] = FULL_VEC_LEN_;
}

// ======================================================================== //
// Layout — choose()
//   Ported from choose() in tpsa.cpp (algorithm unchanged).
//   Returns pointer past the last filled entry.
// ======================================================================== //

static unsigned int* choose_impl(unsigned int n, unsigned int r,
                                  unsigned int* p, unsigned int nv)
{
    if (n == 0 || r == 0) return p;

    std::vector<unsigned int> c(r);
    unsigned int i, j, k;

    for (i = 0; i < r; ++i) c[i] = i;
    c[r - 1] = r - 2;

    while (true) {
        j = r - 1;
        while (j < r && c[j] >= n - r + j && c[j] < n) j--;
        if (j >= r) break;

        ++c[j];
        for (k = j + 1; k < r; ++k) c[k] = c[k - 1] + 1;

        for (size_t iv = 0; iv < r; ++iv)
            ++(*(p + c[iv] - iv));
        p += nv;
    }
    return p;
}

// ======================================================================== //
// Layout — init_base()
//   Ported from init_base() in tpsa.cpp (algorithm unchanged).
//   base_ stores gnv per-row cumulative exponent sums:
//     base_[i*gnv + k] = sum_{l=k}^{gnv-1} exponent_l  (suffix sums).
// ======================================================================== //

void Layout::init_base()
{
    unsigned int nv = gnv_;
    unsigned int nd = gnd_;
    unsigned int fvl = FULL_VEC_LEN_;

    base_ = new unsigned int[nv * fvl];
    std::memset(base_, 0, nv * fvl * sizeof(unsigned int));

    unsigned int* pb = base_;

    // element 0 = constant (all zeros)
    pb += nv;

    for (unsigned int d = 1; d <= nd; ++d) {
        pb = choose_impl(nv + d - 1, d, pb, nv);
    }

    // convert to cumulative (suffix-sum) encoding
    pb = base_;
    for (size_t i = 0; i < fvl; ++i) {
        unsigned int x = 0;
        for (unsigned int k = 0; k < nv; ++k) {
            x += pb[nv - 1 - k];
            pb[nv - 1 - k] = x;
        }
        pb += nv;
    }
}

// ======================================================================== //
// Layout — init_prod_index()
//   Ported from init_prod_index() in tpsa.cpp (algorithm unchanged).
// ======================================================================== //

void Layout::init_prod_index()
{
    unsigned int nv = gnv_;
    unsigned int nd = gnd_;
    unsigned int fvl = FULL_VEC_LEN_;

    // --- Build m and H tables ---
    // m[i][j] and H[i][j] are cumulative binomial-based index helpers.
    std::vector<std::vector<unsigned int>> m(nv + 1, std::vector<unsigned int>(nd + 2, 0));
    H_.assign(nv + 1, std::vector<unsigned int>(nd + 2, 0));

    for (unsigned int i = 0; i < nv + 1; ++i) {
        m[i][0] = H_[i][0] = 0;
        m[i][1] = H_[i][1] = 1;
        for (unsigned int j = 2; j < nd + 2; ++j) {
            m[i][j]   = m[i][j - 1] * (i + j - 2) / (j - 1);
            H_[i][j]  = H_[i][j - 1] * (i + j - 2) / (j - 1);
        }
        for (unsigned int j = 1; j < nd + 2; ++j) {
            m[i][j]  += m[i][j - 1];
            H_[i][j] += H_[i][j - 1];
        }
    }

    // --- Allocate prdidx ---
    // prdidx_ has fvl slots but only indices [1, order_index[nd]) are populated.
    // Indices from order_index[nd] to fvl-1 correspond to highest-order monomials
    // whose products always exceed the truncation order -- they are set to nullptr.
    prdidx_ = new unsigned int*[fvl];
    prdidx_rows_ = fvl;
    for (unsigned int i = 0; i < fvl; ++i) prdidx_[i] = nullptr;

    unsigned int ord = 1;
    const unsigned int* pb = base_;

    unsigned int NS = 0;
    for (size_t i = 1; i < order_index_[nd]; ++i) {
        if (order_index_[ord + 1] <= i) ++ord;
        size_t M = order_index_[nd - ord + 1];
        prdidx_[i] = new unsigned int[M];
        NS += (unsigned int)M;
        for (size_t j = 1; j < M; ++j) {
            prdidx_[i][j] = 0;
            size_t idx = 0;
            for (unsigned int k = 0; k < nv; ++k) {
                idx += m[nv - k][pb[i * nv + k] + pb[j * nv + k]];
            }
            prdidx_[i][j] = (unsigned int)idx;
        }
    }
    tblsize_ = NS;
}

// ======================================================================== //
// Layout — order table (ADOrderTable equivalent)
//   Ported from ADOrderTable::generate_order_table() in tpsa_extend.cc.
// ======================================================================== //

void Layout::generate_order_table()
{
    order_table_.clear();
    order_index_map_.clear();

    unsigned int nv  = gnv_;
    unsigned int fvl = FULL_VEC_LEN_;

    std::vector<int> ord(nv);

    // base_ encodes cumulative suffix sums.
    // The actual exponent of variable j is:
    //   base_[i*nv + j] - base_[i*nv + j + 1]   for j < nv-1
    //   base_[i*nv + nv-1]                        for j == nv-1
    const unsigned int* p = base_;
    for (size_t i = 0; i < fvl; ++i) {
        for (size_t j = 0; j < nv - 1; ++j) {
            ord[j] = (int)(p[j] - p[j + 1]);
        }
        ord[nv - 1] = (int)p[nv - 1];
        p += nv;

        order_table_.push_back(ord);
        order_index_map_.insert({ord, (int)i});
    }

    order_table_valid_ = true;
}

std::vector<int> Layout::orders(unsigned int i) const
{
    if (order_table_valid_) {
        return order_table_.at(i);
    }
    // Fall back: compute from base_ on the fly
    unsigned int nv = gnv_;
    std::vector<int> ord(nv);
    const unsigned int* p = base_ + (size_t)i * nv;
    for (unsigned int j = 0; j < nv - 1; ++j)
        ord[j] = (int)(p[j] - p[j + 1]);
    ord[nv - 1] = (int)p[nv - 1];
    return ord;
}

int Layout::find_index(const std::vector<int>& ord) const
{
    if (!order_table_valid_) {
        throw std::logic_error("Layout::find_index: order table not built — call generate_order_table() first");
    }
    auto it = order_index_map_.find(ord);
    if (it == order_index_map_.end()) {
        throw std::out_of_range("Layout::find_index: exponent vector not found");
    }
    return it->second;
}

// ======================================================================== //
// Layout — change_order / restore_order
//   Ported from ad_change_order / ad_restore_order in tpsa_extend.cc.
// ======================================================================== //

int Layout::change_order(unsigned int new_order)
{
    // Record original order the first time
    if (gnd_record_ == 0) gnd_record_ = gnd_original_;

    if (new_order <= gnd_record_) {
        gnd_ = new_order;
        FULL_VEC_LEN_ = comb_num(gnv_ + gnd_, gnd_);
        return 0;
    }
    return 1;
}

void Layout::restore_order()
{
    if (gnd_record_ > gnd_) {
        gnd_ = gnd_record_;
        FULL_VEC_LEN_ = comb_num(gnv_ + gnd_, gnd_);
    }
}

} // namespace da

// ===========================================================================
// Layout::orders_ref — returns a const reference to the order_table_ entry.
// Builds the order table on demand if not already built.
// ===========================================================================

namespace da {

const std::vector<int>& Layout::orders_ref(unsigned int i)
{
    if (!order_table_valid_) {
        generate_order_table();
    }
    return order_table_.at(i);
}

} // namespace da
