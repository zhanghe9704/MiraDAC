// bind_env.cpp — env-level functions of the current env (plan T1.2) and the
// Env class (plan T7.3).
#include "arith.h"

#include <nanobind/stl/string.h>

#include <functional>
#include <sstream>

namespace {

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

// The env of the last init(), until clear(). C++ has no accessor for it.
da::DAEnv* g_default = nullptr;

// Previous envs of the active `with env:` blocks; they nest, so one stack
// serves every Env.
thread_local std::vector<da::DAEnv*> g_with_stack;

// Env objects hold raw pointers, so no env may be deleted. da_clear() and
// da_destroy_env() delete an env with no live vector and retire one with
// live vectors (keeping its shell until exit); one assigned slot, never
// freed, makes them always retire. The default env's slots held by
// da::base are freed by da_clear()/da_init() before the check, so drop them
// first.
void keep_shell(da::DAEnv* e) {
    if (!e) return;
    if (e == g_default) da::base.base_vecs.clear();
    if (e->live_slots() == 0) e->pool<double>().assign();
}

void check_pool_size(unsigned pool_size) {
    if (pool_size == 0) throw std::invalid_argument("pool_size must be positive");
}

void clear_default() {
    keep_shell(g_default);
    da::da_clear();
    g_default = nullptr;
}

da::DAEnv& live(da::DAEnv& e) {
    if (e.retired()) throw EnvError("DA environment has been cleared");
    return e;
}

template <class DA>
auto import_to(da::DAEnv& dst, const DA& v) {
    EnvGuard g(v.env_);
    live(dst);
    return make_da<DA>([&] { return da::import_to(dst, v); });
}

} // namespace

void bind_env(nb::module_& m) {
    m.def("init", [](unsigned order, unsigned nvars, unsigned pool_size, bool table) {
        check_pool_size(pool_size);
        keep_shell(g_default);
        da::da_init(order, nvars, pool_size, table);
        g_default = &da::da_current_env();
    }, nb::arg("order"), nb::arg("nvars"), nb::arg("pool_size"), nb::arg("table") = false,
       "Create the default DA environment and make it current.\n\n"
       "DA vectors have nvars variables and are truncated at order; the pool\n"
       "holds pool_size vectors. table=True also builds the index/exponent\n"
       "lookup table. A previous default env is cleared first; its remaining\n"
       "vectors raise EnvError on use.");
    m.def("clear", &clear_default,
          "Clear the default env. Its remaining vectors raise EnvError on use.");

    m.def("count", [] { return current_env_or_throw()->pool<double>().count(); },
          "Number of vectors in use in the current env's pool.");
    m.def("remain", [] { return current_env_or_throw()->pool<double>().remain(); },
          "Number of free vectors in the current env's pool.");
    m.def("poolsize", [] { return current_env_or_throw()->pool<double>().poolsize(); },
          "Capacity of the current env's pool.");
    m.def("full_length", [] { return current_env_or_throw()->layout().full_len(); },
          "Number of monomials up to the current order of the current env.");
    m.def("nvars", [] { return current_env_or_throw()->layout().num_vars(); },
          "Number of variables of the current env.");
    m.def("max_order", [] { return original_order(*current_env_or_throw()); },
          "Order the current env was created with.");
    m.def("current_order", [] { return current_env_or_throw()->layout().max_order(); },
          "Order at which the current env truncates now (see change_order()).");

    m.def("get_eps", [] { return da::NDA::eps; },
          "Coefficients below eps in magnitude are dropped (shared by all envs).");
    m.def("set_eps", [](double x) {
        if (!(x > 0)) throw std::invalid_argument("eps must be positive");
        da::NDA::eps = x;
    }, nb::arg("x"), "Set eps (see get_eps()); x must be positive.");

    m.def("change_order", [](unsigned n) {
        return current_env_or_throw()->layout().change_order(n) == 0;
    }, nb::arg("n"),
       "Truncate the current env at order n <= max_order(); False if n is too\n"
       "large. Prefer the order(n) context manager.");
    m.def("restore_order", [] { current_env_or_throw()->layout().restore_order(); },
          "Truncate the current env at its original order (max_order()) again.");

    // Env (T7.3). C++ owns every env and keep_shell() stops any from being
    // deleted, so Python only ever holds references.
    nb::class_<da::DAEnv> e(m, "Env",
        "A DA environment: variables, order and a pool of vectors.\n\n"
        "Env(order, nvars, pool_size, table=False) creates one without changing\n"
        "the current env. `with env:` makes it current for the block; the\n"
        "env= keyword of the constructors creates vectors in it. Every\n"
        "operation on a vector runs in the vector's own env.");
    e.def(nb::new_([](unsigned order, unsigned nvars, unsigned pool_size, bool table) {
        check_pool_size(pool_size);
        struct Restore {
            da::DAEnv* prev = da::da_exchange_env(nullptr);
            ~Restore() { da::da_exchange_env(prev); }
        } r;
        return &da::da_make_env(order, nvars, pool_size, table);
    }), nb::arg("order"), nb::arg("nvars"), nb::arg("pool_size"), nb::arg("table") = false,
        nb::rv_policy::reference, "Create an env; the current env does not change.");
    e.def_static("current", [] {
        da::DAEnv* env = da::da_exchange_env(nullptr);
        da::da_exchange_env(env);
        if (!env) throw EnvError("no current DA environment");
        return env;
    }, nb::rv_policy::reference, "The current env; EnvError if there is none.");
    e.def_static("default", [] {
        if (!g_default) throw EnvError("no default DA environment: call miradac.init() first");
        return g_default;
    }, nb::rv_policy::reference, "The env of the last init(); EnvError after clear().");

    e.def_prop_ro("current_order", [](da::DAEnv& e) { return live(e).layout().max_order(); },
                  "Order at which this env truncates now.");
    e.def_prop_ro("max_order", [](da::DAEnv& e) { return original_order(live(e)); },
                  "Order this env was created with.");
    e.def_prop_ro("nvars", [](da::DAEnv& e) { return live(e).layout().num_vars(); },
                  "Number of variables.");
    e.def_prop_ro("full_length", [](da::DAEnv& e) { return live(e).layout().full_len(); },
                  "Number of monomials up to the current order.");
    e.def_prop_ro("pool_size", [](da::DAEnv& e) { return live(e).pool<double>().poolsize(); },
                  "Capacity of the pool.");
    e.def_prop_ro("count", [](da::DAEnv& e) { return live(e).pool<double>().count(); },
                  "Number of vectors in use.");
    e.def_prop_ro("remain", [](da::DAEnv& e) { return live(e).pool<double>().remain(); },
                  "Number of free vectors.");
    e.def_prop_ro("live_slots", [](da::DAEnv& e) { return live(e).live_slots(); },
                  "Number of vectors in use, numeric and symbolic.");
    e.def_prop_ro("retired", [](const da::DAEnv& e) { return e.retired(); },
                  "True once the env was closed or cleared.");

    e.def("select", [](da::DAEnv& e) { da::da_select_env(live(e)); },
          "Make this env current (until another is selected).");
    e.def("__enter__", [](da::DAEnv& e) -> da::DAEnv& {
        g_with_stack.push_back(da::da_exchange_env(&live(e)));
        return e;
    }, nb::rv_policy::reference, "Make this env current for the `with` block.");
    e.def("__exit__", [](da::DAEnv&, nb::handle, nb::handle, nb::handle) {
        if (g_with_stack.empty()) return;
        da::da_exchange_env(g_with_stack.back());
        g_with_stack.pop_back();
    }, nb::arg("exc_type").none(), nb::arg("exc").none(), nb::arg("tb").none(),
       "Restore the env that was current before the `with` block.");
    e.def("close", [](da::DAEnv& e) {
        if (&e == g_default) return clear_default();
        if (e.retired()) return;
        keep_shell(&e);
        da::da_destroy_env(e);
    }, "Release the pool; the env's vectors raise EnvError on use. Idempotent.");
    e.def("change_order", [](da::DAEnv& e, unsigned n) {
        EnvGuard g(&e);
        return e.layout().change_order(n) == 0;
    }, nb::arg("n"),
       "Truncate this env at order n <= max_order; False if n is too large.\n"
       "Prefer the order(n) context manager.");
    e.def("restore_order", [](da::DAEnv& e) {
        EnvGuard g(&e);
        e.layout().restore_order();
    }, "Truncate this env at its original order (max_order) again.");

    e.def("import_", &import_to<NDA>, nb::arg("v"),
          "Copy v (NDA, CNDA, SDA or CSDA) into this env. ValueError if the\n"
          "two envs have different variables or orders.");
    e.def("import_", [](da::DAEnv& dst, const CNDA& v) {
        EnvGuard g(env_of(v));
        live(dst);
        return make_da<CNDA>([&] {
            return CNDA(da::import_to(dst, da::get_real(v)), da::import_to(dst, da::get_imag(v)));
        });
    }, nb::arg("v"));
#ifdef DA_WITH_SYMBOLIC
    e.def("import_", &import_to<SDA>, nb::arg("v"));
    e.def("import_", [](da::DAEnv& dst, const CSDA& v) {
        EnvGuard g(env_of(v));
        live(dst);
        return make_da<CSDA>([&] {
            return CSDA(da::import_to(dst, da::get_real(v)), da::import_to(dst, da::get_imag(v)));
        });
    }, nb::arg("v"));
    e.def("promote", [](da::DAEnv& dst, const NDA& v) {
        EnvGuard g(v.env_);
        live(dst);
        return make_da<SDA>([&] { return da::promote_to(dst, v); });
    }, nb::arg("v"),
       "promote(v) (NDA to SDA, CNDA to CSDA) with the result in this env.\n"
       "ValueError if the two envs have different variables or orders.");
    e.def("promote", [](da::DAEnv& dst, const CNDA& v) {
        EnvGuard g(env_of(v));
        live(dst);
        return make_da<CSDA>([&] {
            return CSDA(da::promote_to(dst, da::get_real(v)), da::promote_to(dst, da::get_imag(v)));
        });
    }, nb::arg("v"));
#endif

    auto same = [](const da::DAEnv& a, nb::handle b) {
        return nb::isinstance<da::DAEnv>(b) && &nb::cast<const da::DAEnv&>(b) == &a;
    };
    e.def("__eq__", [same](const da::DAEnv& a, nb::handle b) { return same(a, b); },
          nb::arg().none(), "True if value is the same env.");
    e.def("__ne__", [same](const da::DAEnv& a, nb::handle b) { return !same(a, b); },
          nb::arg().none(), "True unless value is the same env.");
    e.def("__hash__", [](const da::DAEnv& e) { return std::hash<const void*>()(&e); },
          "Hash of the env's identity.");
    e.def("__repr__", [](da::DAEnv& e) {
        if (e.retired()) return std::string("Env(retired)");
        std::ostringstream os;
        os << "Env(order=" << original_order(e) << ", nvars=" << e.layout().num_vars()
           << ", pool_size=" << e.pool<double>().poolsize() << ")";
        return os.str();
    }, "Env(order=..., nvars=..., pool_size=...) or Env(retired).");
}
