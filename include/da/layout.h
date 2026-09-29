/**
 * @file layout.h
 * @brief Layout class: type-independent monomial-indexing tables.
 * @details Owns and computes all of the precomputed monomial tables:
 *   - order_index[]  : starting index of each order
 *   - base[]         : exponent vectors for every monomial (cumulative encoding)
 *   - prdidx[][]     : product-index table
 *   - H[][]          : hash table used when building prdidx
 *
 * Also provides the ADOrderTable functionality (index <-> exponent-vector mapping),
 * change_order() / restore_order(), and find_index().
 *
 * This class has no template parameter — it is purely type-independent.
 * It is ported and refactored from:
 *   ref/tpsa/src/tpsa.cpp         (init_order_index, init_prod_index, init_base, ...)
 *   ref/tpsa/src/tpsa_extend.cc   (ADOrderTable, ad_change_order, ad_restore_order)
 */

#pragma once

#include "da/monomial_scheme.h"

#include <vector>
#include <map>
#include <cstdint>
#include <stdexcept>

namespace da {

/**
 * @brief Monomial-indexing tables for a DA environment.
 *
 * A Layout is constructed once from a MonomialScheme and thereafter
 * provides O(1) access to all index-arithmetic needed by the DA kernels.
 *
 * RAII: all heap arrays are owned by Layout and freed in the destructor.
 */
class Layout {
public:
    // ------------------------------------------------------------------ //
    //  Construction / destruction                                          //
    // ------------------------------------------------------------------ //

    /**
     * @brief Build all tables from the given scheme.
     * @param scheme  Describes the monomial set (variables, order).
     * @param build_order_table  If true, also build the order_table for
     *                           index<->exponent-vector lookup (optional
     *                           at construction; can be called later).
     */
    explicit Layout(const MonomialScheme& scheme, bool build_order_table = false);

    ~Layout();

    // Non-copyable (large tables), movable.
    Layout(const Layout&)            = delete;
    Layout& operator=(const Layout&) = delete;
    Layout(Layout&&)                 noexcept;
    Layout& operator=(Layout&&)      noexcept;

    // ------------------------------------------------------------------ //
    //  Basic properties                                                    //
    // ------------------------------------------------------------------ //

    unsigned int num_vars()        const { return gnv_; }
    unsigned int max_order()       const { return gnd_; }
    unsigned int full_len()        const { return FULL_VEC_LEN_; }

    // ------------------------------------------------------------------ //
    //  Raw table accessors (for kernels — resolve to local pointers)       //
    // ------------------------------------------------------------------ //

    /// Starting index of monomials of each order.  Length = gnd+2.
    const unsigned int* order_index() const { return order_index_; }

    /// Flat exponent-vector array (cumulative encoding from tpsa.cpp).
    /// Element i occupies [i*gnv_, (i+1)*gnv_).
    const unsigned int* base() const { return base_; }

    /// prdidx[i][j] = index of the product of monomial i and monomial j.
    /// prdidx[0] == nullptr (constant * anything is a trivial case handled by the kernel).
    /// For i>=1: j in [1, order_index[gnd-ord_i+1]).
    unsigned int** prdidx() const { return prdidx_; }

    /// H[i][j]: cumulative hash table used when computing prdidx.
    const std::vector<std::vector<unsigned int>>& H() const { return H_; }

    // ------------------------------------------------------------------ //
    //  Order-table (index <-> exponent-vector mapping)                     //
    // ------------------------------------------------------------------ //

    /// Build (or rebuild) the order table.
    void generate_order_table();

    /// True if the order table has been built.
    bool has_order_table() const { return order_table_valid_; }

    /**
     * @brief Return the exponent vector of the i-th monomial.
     * @note Requires has_order_table() == true, or falls back to computing
     *       from base[] on the fly.
     */
    std::vector<int> orders(unsigned int i) const;

    /**
     * @brief Return a const reference to the exponent vector of the i-th monomial.
     * @note Requires has_order_table() == true (calls generate_order_table() if not).
     */
    const std::vector<int>& orders_ref(unsigned int i);

    /**
     * @brief Given an exponent vector, return the corresponding index.
     * @note Requires has_order_table() == true.
     * @throws std::out_of_range if the exponent vector is not in the table.
     */
    int find_index(const std::vector<int>& ord) const;

    // ------------------------------------------------------------------ //
    //  Temporary order reduction                                           //
    // ------------------------------------------------------------------ //

    /**
     * @brief Temporarily lower the truncation order.
     * @param new_order  Must be <= original order.
     * @return 0 on success, 1 if new_order > original (no change).
     */
    int  change_order(unsigned int new_order);

    /// Restore the original order given at construction. This is not a stack:
    /// it ignores any lower order set by earlier change_order() calls. To return
    /// to a caller's order, save max_order() and call change_order() with it.
    void restore_order();

    /**
     * @brief Free the monomial tables, keeping the object usable-but-empty.
     *
     * Idempotent. Used by DAEnv::release_memory() to drop the large
     * prdidx/base/order_index allocations when an environment is retired
     * while DAVectors may still reference it.
     */
    void release_tables() noexcept { free_tables(); }

private:
    // ------------------------------------------------------------------ //
    //  Internal helpers                                                    //
    // ------------------------------------------------------------------ //

    static unsigned int comb_num(unsigned int n, unsigned int r);
    static unsigned int gcd_impl(unsigned int a, unsigned int b);

    void init_order_index();
    void init_base();
    void init_prod_index();

    void free_tables();

    // ------------------------------------------------------------------ //
    //  Members (mirror the file-statics in tpsa.cpp / tpsa_extend.cc)    //
    // ------------------------------------------------------------------ //

    unsigned int gnv_;             ///< Number of variables
    unsigned int gnd_;             ///< Current truncation order
    unsigned int gnd_original_;    ///< Original order (before any change_order)
    unsigned int gnd_record_;      ///< Saved order for restore_order (0 = none saved)
    unsigned int FULL_VEC_LEN_;    ///< C(gnv+gnd, gnd)

    unsigned int*  order_index_;   ///< [gnd+2]
    unsigned int*  base_;          ///< [gnv * FULL_VEC_LEN]  (cumulative encoding)
    unsigned int** prdidx_;        ///< [prdidx_rows_], prdidx_[0]==nullptr
    /// Rows allocated in prdidx_. Fixed at construction: FULL_VEC_LEN_ is
    /// lowered by change_order(), so freeing against it would leak the rows
    /// above the reduced order.
    unsigned int   prdidx_rows_;
    unsigned int   tblsize_;       ///< Total entries in prdidx (for diagnostics)

    std::vector<std::vector<unsigned int>> H_;  ///< [gnv+1][gnd+2]

    // Order table (ADOrderTable equivalent)
    bool order_table_valid_;
    std::vector<std::vector<int>>     order_table_;   ///< index -> exponent vector
    std::map<std::vector<int>, int>   order_index_map_; ///< exponent vector -> index
};

} // namespace da
