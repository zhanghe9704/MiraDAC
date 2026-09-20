/**
 * @file env.cpp
 * @brief DAEnv implementation + environment management functions.
 */

#include "da/env.h"

#include <memory>
#include <vector>
#include <stdexcept>

namespace da {

// ======================================================================
// DAEnv construction
// ======================================================================

DAEnv::DAEnv(unsigned order, unsigned num_vars, unsigned pool_size, bool table)
    : layout_(UniformOrder(num_vars, order), table)
{
    unsigned full_len = layout_.full_len();
    pool_d_.reserve(full_len, pool_size);
#ifdef DA_WITH_SYMBOLIC
    pool_e_.reserve(full_len, pool_size);
#endif
}

// ======================================================================
// import<T>: same-type same-layout slot copy
// ======================================================================

template<class T>
unsigned DAEnv::import(const DAEnv& src_env, unsigned src_slot) {
    // Layout compatibility check: full_len and parameter counts must match.
    if (layout_.full_len()  != src_env.layout_.full_len()  ||
        layout_.num_vars()   != src_env.layout_.num_vars()  ||
        layout_.max_order()  != src_env.layout_.max_order()) {
        throw std::invalid_argument(
            "DAEnv::import: source and target layouts are incompatible");
    }
    Pool<T>& dst_pool = pool<T>();
    const Pool<T>& src_pool = src_env.pool<T>();

    unsigned dst_slot = dst_pool.assign();
    unsigned len = src_pool.len(src_slot);
    dst_pool.copy_slot(src_pool.slot(src_slot),
                       dst_pool.slot(dst_slot),
                       len);
    dst_pool.set_len(dst_slot, len);
    return dst_slot;
}

// Explicit instantiations
template unsigned DAEnv::import<double>(const DAEnv&, unsigned);

#ifdef DA_WITH_SYMBOLIC
template unsigned DAEnv::import<SymEngine::Expression>(const DAEnv&, unsigned);
#endif

// ======================================================================
// promote: double -> Expression (same layout)
// ======================================================================

#ifdef DA_WITH_SYMBOLIC
unsigned DAEnv::promote(const DAEnv& src_env, unsigned src_slot) {
    if (layout_.full_len()  != src_env.layout_.full_len()  ||
        layout_.num_vars()   != src_env.layout_.num_vars()  ||
        layout_.max_order()  != src_env.layout_.max_order()) {
        throw std::invalid_argument(
            "DAEnv::promote: source and target layouts are incompatible");
    }
    const Pool<double>& src_pool = src_env.pool<double>();
    Pool<SymEngine::Expression>& dst_pool = pool<SymEngine::Expression>();

    unsigned dst_slot = dst_pool.assign();
    unsigned len = src_pool.len(src_slot);
    const double* src = src_pool.slot(src_slot);
    SymEngine::Expression* dst = dst_pool.slot(dst_slot);
    for (unsigned i = 0; i < len; ++i)
        dst[i] = SymEngine::Expression(src[i]);
    dst_pool.set_len(dst_slot, len);
    return dst_slot;
}
#endif

// ======================================================================
// Thread-local current environment
// ======================================================================

namespace {

// Pointer to the current (thread-local) environment.
thread_local DAEnv* tl_current_env = nullptr;

// Heap-owned environments created by da_make_env / da_init.
// We store them in a vector of unique_ptr so they are properly destroyed.
// da_clear() clears the default env (the last one created by da_init).
static std::unique_ptr<DAEnv> default_env;

} // anonymous namespace

DAEnv& da_current_env() {
    if (!tl_current_env) {
        throw std::runtime_error(
            "da::da_current_env: no DA environment has been selected. "
            "Call da_init() first.");
    }
    return *tl_current_env;
}

void da_select_env(DAEnv& env) {
    tl_current_env = &env;
}

DAEnv& da_make_env(unsigned order, unsigned num_vars,
                   unsigned pool_size, bool table) {
    // Create on the heap; the caller is responsible for keeping it alive.
    // For the simple multi-env case, we return a reference to a static
    // vector-managed object.  The user must keep their own handle or use
    // the returned reference while the env is in scope.
    // We allocate with new and return a reference; ownership stays with
    // the caller (or they can let it go and use da_clear later).
    // Note: this is intentionally minimal per the design plan (Stage 3).
    auto* env = new DAEnv(order, num_vars, pool_size, table);
    tl_current_env = env;
    return *env;
}

// ======================================================================
// Classic single-environment entry points
// ======================================================================
// da_init() and da_clear() are defined in src/base.cpp so they can
// call base.set_base() without a circular include dependency.

int da_init_env_only(unsigned order, unsigned num_vars,
                     unsigned pool_size, bool table) {
    default_env = std::make_unique<DAEnv>(order, num_vars, pool_size, table);
    tl_current_env = default_env.get();
    return 0;
}

void da_clear_env_only() {
    if (default_env) {
        if (tl_current_env == default_env.get()) {
            tl_current_env = nullptr;
        }
        default_env.reset();
    }
}

} // namespace da
