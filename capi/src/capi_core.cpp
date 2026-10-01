// capi_core.cpp — versions, error handling and env functions of the C API
// (plan T0.1, T1.2).
#include "common.h"

#include <string>

using namespace mdac;

namespace {

thread_local std::string g_last_error;

mdac_status fail(mdac_status s, const char* what) noexcept {
    try {
        g_last_error = what;
    } catch (...) {
    }
    return s;
}

// The order given at construction. Layout keeps it private, but the pool
// slot length is fixed at construction to C(nvars + order, order).
unsigned original_order(da::DAEnv& e) {
    const da::Layout& l = e.layout();
    const unsigned full = e.pool<double>().full_len();
    unsigned n = l.max_order();
    for (unsigned len = l.full_len(); len < full; ++n)
        len = len * (l.num_vars() + n + 1) / (n + 1);
    return n;
}

// The env of the last mdac_init(), until mdac_clear(). C++ has no accessor.
da::DAEnv* g_default = nullptr;

// Env handles are raw pointers, so no env may be deleted. da_clear() and
// da_destroy_env() delete an env with no live vector and retire one with
// live vectors (keeping its shell until exit); one assigned slot, never
// freed, makes them always retire. The default env's slots held by da::base
// are freed by da_clear()/da_init() before the check, so drop them first.
void keep_shell(da::DAEnv* e) {
    if (!e) return;
    if (e == g_default) da::base.base_vecs.clear();
    if (e->live_slots() == 0) e->pool<double>().assign();
}

void check_pool_size(unsigned poolsize) {
    if (poolsize == 0) throw std::invalid_argument("poolsize must be positive");
}

void clear_default() {
    keep_shell(g_default);
    da::da_clear();
    g_default = nullptr;
}

da::DAEnv& live(mdac_env* h) {
    da::DAEnv* e = env(h);
    if (!e) throw EnvError("no DA environment");
    if (e->retired()) throw EnvError("DA environment has been cleared");
    return *e;
}

} // namespace

mdac_status mdac::current_exception_status() noexcept {
    try {
        throw;
    } catch (const da::PoolExhausted& e) {
        return fail(MDAC_ERR_POOL, e.what());
    } catch (const EnvError& e) {
        return fail(MDAC_ERR_ENV, e.what());
    } catch (const std::invalid_argument& e) {
        return fail(MDAC_ERR_VALUE, e.what());
    } catch (const std::domain_error& e) {
        return fail(MDAC_ERR_VALUE, e.what());
    } catch (const std::out_of_range& e) {
        return fail(MDAC_ERR_INDEX, e.what());
    } catch (const std::logic_error& e) {  // da::check_env
        return fail(MDAC_ERR_ENV, e.what());
    } catch (const std::exception& e) {
        return fail(MDAC_ERR_RUNTIME, e.what());
    } catch (...) {
        return fail(MDAC_ERR_RUNTIME, "unknown C++ exception");
    }
}

extern "C" {

int mdac_abi_version(void) { return MDAC_ABI_VERSION; }

const char* mdac_version(void) { return MDAC_VERSION; }

const char* mdac_last_error(void) { return g_last_error.c_str(); }

mdac_status mdac_init(unsigned order, unsigned nvars, unsigned poolsize, int table) {
    MDAC_TRY {
        check_pool_size(poolsize);
        keep_shell(g_default);
        da::da_init(order, nvars, poolsize, table != 0);
        g_default = &da::da_current_env();
    } MDAC_CATCH
}

mdac_status mdac_clear(void) {
    MDAC_TRY { clear_default(); } MDAC_CATCH
}

mdac_status mdac_env_current(mdac_env** out) {
    MDAC_TRY {
        da::DAEnv* e = da::da_exchange_env(nullptr);
        da::da_exchange_env(e);
        if (!e) throw EnvError("no current DA environment");
        *out = handle(e);
    } MDAC_CATCH
}

mdac_status mdac_env_default(mdac_env** out) {
    MDAC_TRY {
        if (!g_default) throw EnvError("no default DA environment: call mdac_init() first");
        *out = handle(g_default);
    } MDAC_CATCH
}

mdac_status mdac_env_make(unsigned order, unsigned nvars, unsigned poolsize, int table,
                          mdac_env** out) {
    MDAC_TRY {
        check_pool_size(poolsize);
        struct Restore {
            da::DAEnv* prev = da::da_exchange_env(nullptr);
            ~Restore() { da::da_exchange_env(prev); }
        } r;
        *out = handle(&da::da_make_env(order, nvars, poolsize, table != 0));
    } MDAC_CATCH
}

mdac_status mdac_env_select(mdac_env* e) {
    MDAC_TRY { da::da_select_env(live(e)); } MDAC_CATCH
}

mdac_status mdac_env_exchange(mdac_env* e, mdac_env** prev) {
    MDAC_TRY { *prev = handle(da::da_exchange_env(env(e))); } MDAC_CATCH
}

mdac_status mdac_env_close(mdac_env* h) {
    MDAC_TRY {
        da::DAEnv* e = env(h);
        if (!e) throw EnvError("no DA environment");
        if (e == g_default) {
            clear_default();
        } else if (!e->retired()) {
            keep_shell(e);
            da::da_destroy_env(*e);
        }
    } MDAC_CATCH
}

#define MDAC_ENV_QUERY(name, expr)                                    \
    mdac_status mdac_env_##name(mdac_env* h, unsigned* out) {         \
        MDAC_TRY { da::DAEnv& e = live(h); *out = (expr); } MDAC_CATCH \
    }

MDAC_ENV_QUERY(order, e.layout().max_order())
MDAC_ENV_QUERY(max_order, original_order(e))
MDAC_ENV_QUERY(nvars, e.layout().num_vars())
MDAC_ENV_QUERY(full_length, e.layout().full_len())
MDAC_ENV_QUERY(poolsize, e.pool<double>().poolsize())
MDAC_ENV_QUERY(count, e.pool<double>().count())
MDAC_ENV_QUERY(remain, e.pool<double>().remain())

int mdac_env_retired(const mdac_env* e) {
    return reinterpret_cast<const da::DAEnv*>(e)->retired() ? 1 : 0;
}

mdac_status mdac_env_change_order(mdac_env* h, unsigned n, int* ok) {
    MDAC_TRY {
        EnvGuard g(env(h));
        *ok = env(h)->layout().change_order(n) == 0;
    } MDAC_CATCH
}

mdac_status mdac_env_restore_order(mdac_env* h) {
    MDAC_TRY {
        EnvGuard g(env(h));
        env(h)->layout().restore_order();
    } MDAC_CATCH
}

double mdac_get_eps(void) { return da::NDA::eps; }

mdac_status mdac_set_eps(double x) {
    MDAC_TRY {
        if (!(x > 0)) throw std::invalid_argument("eps must be positive");
        da::NDA::eps = x;
    } MDAC_CATCH
}

}
