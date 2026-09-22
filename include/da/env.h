/**
 * @file env.h
 * @brief DAEnv: container for a Layout + Pool(s), and environment management.
 *
 * @details A DAEnv encapsulates all formerly-global state:
 *   - one Layout (type-independent monomial tables), and
 *   - one Pool<double>
 *   - one Pool<SymEngine::Expression>  [only when DA_WITH_SYMBOLIC]
 *
 * The classic single-environment API (da_init / da_clear) manages a default
 * environment selected via a thread-local current-env pointer.
 */

#pragma once

#include "da/layout.h"
#include "da/pool.h"
#include "da/monomial_scheme.h"

#include <memory>
#include <stdexcept>

#ifdef DA_WITH_SYMBOLIC
#  include <symengine/expression.h>
#endif

namespace da {

// ======================================================================
// DAEnv
// ======================================================================

/**
 * @brief Owns a Layout and the coefficient pool(s) for one DA environment.
 *
 * Non-copyable (large tables); movable.
 */
class DAEnv {
public:
    // ------------------------------------------------------------------ //
    //  Construction                                                        //
    // ------------------------------------------------------------------ //

    /**
     * @brief Build a full DA environment.
     *
     * @param order      Truncation order.
     * @param num_vars   Number of independent variables.
     * @param pool_size  Number of slots to pre-allocate in each pool.
     * @param table      If true, generate the order table (index<->exponent map).
     */
    DAEnv(unsigned order, unsigned num_vars, unsigned pool_size,
          bool table = false);

    ~DAEnv() = default;

    // Non-copyable, movable
    DAEnv(const DAEnv&)            = delete;
    DAEnv& operator=(const DAEnv&) = delete;
    DAEnv(DAEnv&&)                 = default;
    DAEnv& operator=(DAEnv&&)      = default;

    // ------------------------------------------------------------------ //
    //  Accessors                                                           //
    // ------------------------------------------------------------------ //

    Layout&       layout()       { return layout_; }
    const Layout& layout() const { return layout_; }

    // ------------------------------------------------------------------ //
    //  Pool accessors — specialised via if constexpr below                 //
    // ------------------------------------------------------------------ //

    /**
     * @brief Return the pool for coefficient type T.
     *
     * Specializations exist for double (and SymEngine::Expression when enabled).
     * Any other T produces a compile-time error.
     */
    template<class T>
    Pool<T>& pool();

    template<class T>
    const Pool<T>& pool() const;

    // ------------------------------------------------------------------ //
    //  Import / convert                                                    //
    // ------------------------------------------------------------------ //

    /**
     * @brief Copy a slot from another env (or this env) of the same type.
     *
     * Same-layout transfer: allocates a new slot in THIS env's pool,
     * copies the data, and returns the new slot index.
     *
     * @throws std::invalid_argument if layouts differ.
     */
    template<class T>
    unsigned import(const DAEnv& src_env, unsigned src_slot);

    // ------------------------------------------------------------------ //
    //  Retirement                                                          //
    // ------------------------------------------------------------------ //

    /**
     * @brief Number of pool slots still held, across all coefficient types.
     *
     * At teardown time (no kernel running) this is the number of live
     * DAVectors that still point at this environment.
     */
    unsigned live_slots() const noexcept;

    /**
     * @brief Release the pools and monomial tables, keeping the object alive.
     *
     * A DAVector must not outlive its DAEnv: ~DAVector reads poolsize()
     * through env_, so destroying the environment first is a use-after-free.
     * The reference library kept its pool in never-destroyed globals, so
     * clearing before a vector went out of scope was harmless there.
     *
     * Retiring restores that property. Everything expensive is freed here,
     * while SymEngine is still initialized, and the empty shell stays valid
     * so late destructors read poolsize() == 0 and skip their free().
     */
    void release_memory() noexcept;

    /// True once release_memory() has run.
    bool retired() const noexcept { return retired_; }

#ifdef DA_WITH_SYMBOLIC
    /**
     * @brief Promote a double slot to a SymEngine::Expression slot.
     *
     * Same-layout only; each element is promoted via
     * Expression(double_value).
     *
     * @throws std::invalid_argument if layouts differ.
     */
    unsigned promote(const DAEnv& src_env, unsigned src_slot);
#endif

private:
    Layout        layout_;
    Pool<double>  pool_d_;
#ifdef DA_WITH_SYMBOLIC
    Pool<SymEngine::Expression> pool_e_;
#endif
    bool          retired_ = false;
};

// ---------------------------------------------------------------------- //
// Pool<T> accessor specializations (inline, in header)                    //
// ---------------------------------------------------------------------- //

template<>
inline Pool<double>& DAEnv::pool<double>() { return pool_d_; }

template<>
inline const Pool<double>& DAEnv::pool<double>() const { return pool_d_; }

#ifdef DA_WITH_SYMBOLIC
template<>
inline Pool<SymEngine::Expression>& DAEnv::pool<SymEngine::Expression>() { return pool_e_; }

template<>
inline const Pool<SymEngine::Expression>& DAEnv::pool<SymEngine::Expression>() const { return pool_e_; }
#endif

// ======================================================================
// Environment management (thread-local current env)
// ======================================================================

/**
 * @brief Return the current (thread-local) DA environment.
 * @throws std::runtime_error if no environment has been selected.
 */
DAEnv& da_current_env();

/**
 * @brief Select env as the current (thread-local) DA environment.
 * @param env  Must outlive any DAVectors created while it is current.
 */
void da_select_env(DAEnv& env);

/**
 * @brief Create a new heap-allocated environment and select it as current.
 *
 * The returned reference is owned by an internal static list.
 * Lifetime: until da_clear() is called for it.
 */
DAEnv& da_make_env(unsigned order, unsigned num_vars,
                   unsigned pool_size, bool table = false);

/**
 * @brief Destroy an environment returned by da_make_env().
 *
 * Prefer this over `delete &env`. If no DAVector still references env, it
 * is deleted outright. If some do, env is retired instead (see
 * DAEnv::release_memory()): its memory is freed now and the empty shell is
 * kept until process exit, so those vectors' destructors stay valid.
 */
void da_destroy_env(DAEnv& env);

// ======================================================================
// Classic single-environment entry points
// ======================================================================

/**
 * @brief Initialize (or re-initialize) the default DA environment.
 *
 * Creates a new DAEnv with the given parameters, selects it as current.
 * Reference default: table = false (numeric), table = true (symbolic).
 * Here we follow the numeric default (false).
 *
 * @return 0 on success (mirrors the reference return type).
 */
int da_init(unsigned order, unsigned num_vars,
            unsigned pool_size, bool table = false);

/**
 * @brief Tear down the default environment created by da_init().
 *
 * Safe to call even if da_init() has not been called (no-op).
 * After this call, da_current_env() will throw until da_init() is
 * called again.
 */
void da_clear();

// Internal helpers called by da_init/da_clear in base.cpp
// (These set up env state without touching the bases object)
int  da_init_env_only(unsigned order, unsigned num_vars,
                      unsigned pool_size, bool table = false);
void da_clear_env_only();

// ======================================================================
// Environment guard helpers
// ======================================================================

/**
 * @brief Return true if both pointers point to the same DAEnv.
 *
 * Intended for use by DAVector binary operators.
 */
inline bool same_env(const DAEnv* a, const DAEnv* b) noexcept {
    return a == b;
}

/**
 * @brief Throw std::logic_error if a and b are different environments.
 *
 * The check is compiled-in only when DA_CHECK_ENV != 0.
 * When DA_CHECK_ENV == 0 the function is a no-op.
 */
// Default the guard ON when nothing defined it. The build system passes
// DA_CHECK_ENV as a PUBLIC definition, but code compiled with a bare
// -Iinclude would otherwise get DA_CHECK_ENV == 0 from the preprocessor's
// "an undefined identifier evaluates to 0" rule and silently lose the check.
#ifndef DA_CHECK_ENV
#define DA_CHECK_ENV 1
#endif

inline void check_env(const DAEnv* a, const DAEnv* b) {
#if DA_CHECK_ENV
    if (a != b) {
        throw std::logic_error(
            "da::check_env: operands belong to different environments");
    }
#else
    (void)a; (void)b;
#endif
}

} // namespace da
