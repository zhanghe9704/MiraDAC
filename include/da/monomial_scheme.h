/**
 * @file monomial_scheme.h
 * @brief MonomialScheme abstract interface and UniformOrder implementation.
 * @details Describes the set of admissible exponent vectors (monomials) in the DA.
 *          The layout layer uses a MonomialScheme to compute the monomial tables.
 */

#pragma once

#include <cstdint>

namespace da {

/**
 * @brief Abstract interface describing the admissible monomial set.
 *
 * A MonomialScheme determines which exponent vectors are present in the
 * truncated power series, and in what order.  It is consulted only at
 * table-build time (inside the Layout constructor); the kernels never
 * call into it at runtime.
 */
class MonomialScheme {
public:
    virtual ~MonomialScheme() = default;

    /// Number of independent variables.
    virtual unsigned int num_vars() const = 0;

    /// Maximum total order (i.e. the truncation order).
    virtual unsigned int max_total_order() const = 0;

    /**
     * @brief Total number of monomials up to and including max_total_order().
     *
     * For uniform order this equals C(nv + nd, nd).
     */
    virtual unsigned int full_len() const = 0;
};

/**
 * @brief Uniform-order monomial scheme: all monomials with total degree <= nd.
 *
 * This is the only MonomialScheme implementation in v1.
 * full_len() == C(nv + nd, nd).
 */
class UniformOrder : public MonomialScheme {
public:
    /**
     * @param nv  Number of independent variables.
     * @param nd  Maximum total order (truncation order).
     */
    UniformOrder(unsigned int nv, unsigned int nd);

    unsigned int num_vars()        const override { return nv_; }
    unsigned int max_total_order() const override { return nd_; }
    unsigned int full_len()        const override;

private:
    unsigned int nv_;
    unsigned int nd_;
};

} // namespace da
