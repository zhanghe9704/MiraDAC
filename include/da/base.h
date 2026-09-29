/**
 * @file base.h
 * @brief Base struct and da::base accessor.
 *
 * @details In the reference, the bases global object is named "da" (da[i]).
 *   Since our namespace is also "da", the bases object is renamed to "base"
 *   per plan A.3: accessed as da::base[i].
 *
 * The Base struct holds std::vector<NDA> and provides operator[].
 * The global "base" object is bound to the current env at the first call to
 * set_base() (called by da_init()).
 */

#pragma once

#include "da/davector.h"

namespace da {

/**
 * @brief Holds the base DA vectors (one per variable).
 *
 * The i-th base is a DA vector with constant = 0 and the i-th linear
 * coefficient = 1.  Access via da::base[i].
 */
struct Base {
    std::vector<NDA> base_vecs;

    Base() {}
    explicit Base(unsigned int n) { set_base(n); }

    /// Initialize n base vectors using ad_var.
    void set_base(unsigned int n);

    /// Initialize using the current env's dimension.
    void set_base();

    /// Access the i-th base.
    const NDA& operator[](unsigned int i) const { return base_vecs.at(i); }
    NDA& operator[](unsigned int i) { return base_vecs.at(i); }
};

/// Global bases accessor — da::base[i] gives the i-th DA base.
extern Base base;

/**
 * @brief Build the i-th base vector in the current env.
 *
 * Unlike da::base[i], which belongs to the default env created by da_init(),
 * this works in whichever env is current. Library code must use it instead of
 * da::base so that it stays correct in envs made by da_make_env().
 */
NDA da_base(unsigned int i);

} // namespace da
