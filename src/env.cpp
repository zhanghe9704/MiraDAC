/**
 * @file env.cpp
 * @brief DAEnv implementation + environment management functions.
 */

#include "da/env.h"

#include <memory>
#include <vector>
#include <stdexcept>
#include <utility>   // std::move

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
// Retirement
// ======================================================================

unsigned DAEnv::live_slots() const noexcept {
    unsigned n = pool_d_.count();
#ifdef DA_WITH_SYMBOLIC
    n += pool_e_.count();
#endif
    return n;
}

void DAEnv::release_memory() noexcept {
    // Release the Expression pool while SymEngine is still initialized.
    // Doing it here, rather than leaving it to static destruction, is what
    // keeps the retired shells inert at process exit.
#ifdef DA_WITH_SYMBOLIC
    pool_e_.release();
#endif
    pool_d_.release();          // sets poolsize() to 0; free() becomes a no-op
    layout_.release_tables();   // the large prdidx / base / order_index blocks
    retired_ = true;
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
    // Same rule as promote(const NDA&): an integer value becomes an exact integer.
    for (unsigned i = 0; i < len; ++i) {
        long iv = static_cast<long>(src[i]);
        if (src[i] == static_cast<double>(iv))
            dst[i] = SymEngine::Expression(iv);
        else
            dst[i] = SymEngine::Expression(src[i]);
    }
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

// Environments cleared while DAVectors still referenced them.
// release_memory() has already handed back everything they held, so each
// entry is an empty shell (a few hundred bytes) kept valid until process
// exit. Owning them here rather than leaking them keeps ASan runs clean.
std::vector<std::unique_ptr<DAEnv>>& retired_envs() {
    static std::vector<std::unique_ptr<DAEnv>> envs;
    return envs;
}

void retire(std::unique_ptr<DAEnv> env) {
    env->release_memory();
    retired_envs().push_back(std::move(env));
}

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

DAEnv* da_exchange_env(DAEnv* env) noexcept {
    DAEnv* prev = tl_current_env;
    tl_current_env = env;
    return prev;
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
    // Re-initializing drops the previous default env, so it goes through the
    // same retirement path as da_clear(); otherwise DAVectors still pointing
    // at the old env would be left dangling.
    da_clear_env_only();
    default_env = std::make_unique<DAEnv>(order, num_vars, pool_size, table);
    tl_current_env = default_env.get();
    return 0;
}

void da_clear_env_only() {
    if (!default_env) return;
    if (tl_current_env == default_env.get()) {
        tl_current_env = nullptr;
    }
    if (default_env->live_slots() == 0) {
        default_env.reset();              // nothing points at it: free it now
    } else {
        // DAVectors are still alive. Free this env's memory but keep the
        // shell valid so their destructors can read poolsize() == 0.
        retire(std::move(default_env));
    }
}

void da_destroy_env(DAEnv& env) {
    if (env.live_slots() == 0) {
        delete &env;
    } else {
        env.release_memory();
        retired_envs().emplace_back(&env);
    }
}

} // namespace da
