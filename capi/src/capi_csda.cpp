// capi_csda.cpp — CSDA (std::complex<SDA>), CSDA lists, promote, evaluate and
// cd_composition for T = Expression in the C API (plan T6.1). Without
// DA_WITH_SYMBOLIC every function is a stub (MDAC_SYM, MDAC_SYM_ELSE).
#include "common.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifdef DA_WITH_SYMBOLIC
using namespace mdac;
using da::get_imag;
using da::get_real;
using Complex = std::complex<double>;

namespace {

da::DAEnv* env_of(const NDA& v) { return v.env_; }
da::DAEnv* env_of(const SDA& v) { return v.env_; }
da::DAEnv* env_of(const CNDA& v) { return get_real(v).env_; }
da::DAEnv* env_of(const CSDA& v) { return get_real(v).env_; }

// The env shared by every DA argument.
template <class A, class... Rest>
da::DAEnv* common_env(const A& a, const Rest&... rest) {
    da::DAEnv* e = env_of(a);
    if (((env_of(rest) != e) || ...)) throw EnvError("DA vectors belong to different environments");
    return e;
}

// out = t part by part, keeping out's slots (a move assignment would swap them).
void assign(CSDA& out, const CSDA& t) {
    get_real(out) = get_real(t);
    get_imag(out) = get_imag(t);
}

// A new R f() in the env of the DA arguments ds.
template <class R, class H, class F, class... D>
mdac_status make(H** out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ds...));
        *out = mdac::handle(new R(f()));
    } MDAC_CATCH
}

// out = f() in out's own slots; out may be one of ds.
template <class F, class... D>
mdac_status into(mdac_csda* out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ref(out), ds...));
        assign(ref(out), f());
    } MDAC_CATCH
}

template <class L>
void same_env(da::DAEnv* e, const L& l) {
    for (const auto& v : l)
        if (env_of(v) != e) throw EnvError("DA vectors belong to different environments");
}

// cd_composition(m, v, out) as a new list; C++ only asserts the checks.
template <class M, class V>
mdac_status compose(const M& m, const V& v, mdac_csdalist** out) {
    MDAC_TRY {
        if (m.empty()) { *out = handle(new CSDAList()); return MDAC_OK; }
        da::DAEnv* e = env_of(m[0]);
        EnvGuard g(e);
        if (v.size() != e->layout().num_vars())
            throw std::invalid_argument("compose: the arguments must be nvars vectors");
        same_env(e, m);
        same_env(e, v);
        auto o = std::make_unique<CSDAList>(m.size());
        da::cd_composition(const_cast<M&>(m), const_cast<V&>(v), *o);
        *out = handle(o.release());
    } MDAC_CATCH
}

} // namespace
#endif // DA_WITH_SYMBOLIC

extern "C" {

// ---- Lifecycle and inspection -------------------------------------------------

mdac_status mdac_csda_new(const mdac_sda* re, const mdac_sda* im, mdac_csda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(im ? common_env(ref(re), ref(im)) : env_of(ref(re)));
        *out = handle(new CSDA(ref(re), im ? ref(im) : SDA()));
    } MDAC_CATCH)
}

mdac_status mdac_csda_new_z(mdac_env* e, double re, double im, mdac_csda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env(e));
        *out = handle(new CSDA(SDA(re), SDA(im)));
    } MDAC_CATCH)
}

mdac_status mdac_csda_promote(const mdac_cnda* v, mdac_csda** out) {
    MDAC_SYM(return make<CSDA>(out, [&] { return da::promote(ref(v)); }, ref(v));)
}

mdac_status mdac_csda_copy(const mdac_csda* v, mdac_csda** out) {
    MDAC_SYM(return make<CSDA>(out, [&] { return ref(v); }, ref(v));)
}

void mdac_csda_free(mdac_csda* v) { MDAC_SYM_ELSE((void)v;, delete &ref(v);) }

mdac_env* mdac_csda_env(const mdac_csda* v) {
    MDAC_SYM_ELSE((void)v; return nullptr;, return handle(env_of(ref(v)));)
}

mdac_status mdac_csda_import(mdac_env* e, const mdac_csda* v, mdac_csda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env_of(ref(v)));
        da::DAEnv& d = live(e);
        *out = handle(new CSDA(da::import_to(d, get_real(ref(v))), da::import_to(d, get_imag(ref(v)))));
    } MDAC_CATCH)
}

mdac_status mdac_csda_real(const mdac_csda* v, mdac_sda** out) {
    MDAC_SYM(return make<SDA>(out, [&] { return get_real(ref(v)); }, ref(v));)
}

mdac_status mdac_csda_imag(const mdac_csda* v, mdac_sda** out) {
    MDAC_SYM(return make<SDA>(out, [&] { return get_imag(ref(v)); }, ref(v));)
}

mdac_status mdac_csda_set_real(mdac_csda* v, const mdac_sda* x) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(common_env(ref(v), ref(x)));
        get_real(ref(v)) = ref(x);
    } MDAC_CATCH)
}

mdac_status mdac_csda_set_imag(mdac_csda* v, const mdac_sda* x) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(common_env(ref(v), ref(x)));
        get_imag(ref(v)) = ref(x);
    } MDAC_CATCH)
}

mdac_status mdac_csda_to_string(const mdac_csda* v, char* buf, size_t cap, size_t* n) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env_of(ref(v)));
        std::ostringstream os;
        os << ref(v);
        const std::string s = os.str();
        *n = s.size() + 1;
        if (cap > 0) {
            const size_t m = std::min(cap - 1, s.size());
            std::memcpy(buf, s.data(), m);
            buf[m] = '\0';
        }
    } MDAC_CATCH)
}

// ---- Arithmetic ------------------------------------------------------------------
// One allocating and one _into form per operand pair: the arguments (after
// the operands' declarations), the C++ expression, and the DA operands.

#define MDAC_CSDA_PAIR(name, expr, ds, ...)                                                  \
    mdac_status name(__VA_ARGS__, mdac_csda** out) {                                         \
        MDAC_SYM(return make<CSDA>(out, [&] { return expr; }, MDAC_CSDA_UNPAREN ds);)        \
    }                                                                                        \
    mdac_status name##_into(mdac_csda* out, __VA_ARGS__) {                                   \
        MDAC_SYM(return into(out, [&] { return expr; }, MDAC_CSDA_UNPAREN ds);)              \
    }

#define MDAC_CSDA_UNPAREN(...) __VA_ARGS__

#define MDAC_CSDA_BINOP(name, op)                                                            \
    MDAC_CSDA_PAIR(mdac_csda_##name, ref(a) op ref(b), (ref(a), ref(b)), const mdac_csda* a, const mdac_csda* b) \
    MDAC_CSDA_PAIR(mdac_csda_##name##_s, ref(a) op ref(b), (ref(a), ref(b)), const mdac_csda* a, const mdac_sda* b) \
    MDAC_CSDA_PAIR(mdac_csda_s##name, ref(a) op ref(b), (ref(a), ref(b)), const mdac_sda* a, const mdac_csda* b) \
    MDAC_CSDA_PAIR(mdac_csda_##name##_e, ref(a) op ref(x), (ref(a)), const mdac_csda* a, const mdac_expr* x) \
    MDAC_CSDA_PAIR(mdac_csda_e##name, ref(x) op ref(a), (ref(a)), const mdac_expr* x, const mdac_csda* a) \
    MDAC_CSDA_PAIR(mdac_csda_##name##_d, ref(a) op x, (ref(a)), const mdac_csda* a, double x) \
    MDAC_CSDA_PAIR(mdac_csda_d##name, x op ref(a), (ref(a)), double x, const mdac_csda* a)  \
    MDAC_CSDA_PAIR(mdac_csda_##name##_z, ref(a) op Complex(re, im), (ref(a)), const mdac_csda* a, double re, double im) \
    MDAC_CSDA_PAIR(mdac_csda_z##name, Complex(re, im) op ref(a), (ref(a)), double re, double im, const mdac_csda* a) \
    MDAC_CSDA_PAIR(mdac_sda_##name##_z, ref(a) op Complex(re, im), (ref(a)), const mdac_sda* a, double re, double im) \
    MDAC_CSDA_PAIR(mdac_sda_z##name, Complex(re, im) op ref(a), (ref(a)), double re, double im, const mdac_sda* a)

MDAC_CSDA_BINOP(add, +)
MDAC_CSDA_BINOP(sub, -)
MDAC_CSDA_BINOP(mul, *)
MDAC_CSDA_BINOP(div, /)

MDAC_CSDA_PAIR(mdac_csda_pow_i, da::pow(ref(a), n), (ref(a)), const mdac_csda* a, int n)
MDAC_CSDA_PAIR(mdac_csda_pow_d, da::pow(ref(a), x), (ref(a)), const mdac_csda* a, double x)

mdac_status mdac_csda_neg(const mdac_csda* a, mdac_csda** out) {
    MDAC_SYM(return make<CSDA>(out, [&] { return -ref(a); }, ref(a));)
}

mdac_status mdac_csda_neg_into(mdac_csda* out, const mdac_csda* a) {
    MDAC_SYM(return into(out, [&] { return -ref(a); }, ref(a));)
}

// ---- Math functions --------------------------------------------------------------

#define MDAC_CSDA_FUNC(f)                                                                    \
    mdac_status mdac_csda_##f(const mdac_csda* a, mdac_csda** out) {                        \
        MDAC_SYM(return make<CSDA>(out, [&] { return da::f(ref(a)); }, ref(a));)             \
    }                                                                                        \
    mdac_status mdac_csda_##f##_into(mdac_csda* out, const mdac_csda* a) {                  \
        MDAC_SYM(return into(out, [&] { return da::f(ref(a)); }, ref(a));)                   \
    }

MDAC_CSDA_FUNC(sqrt)
MDAC_CSDA_FUNC(exp)
MDAC_CSDA_FUNC(log)
MDAC_CSDA_FUNC(asin)
MDAC_CSDA_FUNC(acos)
MDAC_CSDA_FUNC(atan)
MDAC_CSDA_FUNC(asinh)
MDAC_CSDA_FUNC(acosh)
MDAC_CSDA_FUNC(atanh)

mdac_status mdac_csda_abs(const mdac_csda* a, mdac_sda** out) {
    MDAC_SYM(return make<SDA>(out, [&] { return da::abs(ref(a)); }, ref(a));)
}

// The vector overload of da::evaluate (interop.h), as mdac_sda_evaluate.
mdac_status mdac_csda_evaluate(const mdac_csda* v, size_t n, const mdac_expr* const* keys,
                               const double* vals, mdac_cnda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env_of(ref(v)));
        std::vector<SymEngine::RCP<const SymEngine::Basic>> syms;
        for (size_t i = 0; i < n; ++i) syms.push_back(ref(keys[i]).get_basic());
        const std::vector<double> x(vals, vals + n);
        *out = handle(new CNDA(da::evaluate(ref(v), syms, x)));
    } MDAC_CATCH)
}

// ---- Lists -------------------------------------------------------------------------

mdac_status mdac_csdalist_new(mdac_csdalist** out) {
    MDAC_SYM(MDAC_TRY { *out = handle(new CSDAList()); } MDAC_CATCH)
}

mdac_status mdac_csdalist_from(const mdac_csda* const* vs, size_t n, mdac_csdalist** out) {
    MDAC_SYM(MDAC_TRY {
        auto l = std::make_unique<CSDAList>();
        if (n > 0) {
            EnvGuard g(env_of(ref(vs[0])));
            for (size_t i = 1; i < n; ++i) common_env(ref(vs[0]), ref(vs[i]));
            l->reserve(n);
            for (size_t i = 0; i < n; ++i) l->push_back(ref(vs[i]));
        }
        *out = handle(l.release());
    } MDAC_CATCH)
}

void mdac_csdalist_free(mdac_csdalist* l) { MDAC_SYM_ELSE((void)l;, delete &ref(l);) }

size_t mdac_csdalist_length(const mdac_csdalist* l) {
    MDAC_SYM_ELSE((void)l; return 0;, return ref(l).size();)
}

mdac_env* mdac_csdalist_env(const mdac_csdalist* l) {
    MDAC_SYM_ELSE((void)l; return nullptr;,
                  return ref(l).empty() ? nullptr : handle(env_of(ref(l)[0]));)
}

mdac_status mdac_csdalist_get(const mdac_csdalist* l_, size_t i, mdac_csda** out) {
    MDAC_SYM(MDAC_TRY {
        const CSDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(env_of(l[i]));
        *out = handle(new CSDA(l[i]));
    } MDAC_CATCH)
}

mdac_status mdac_csdalist_set(mdac_csdalist* l_, size_t i, const mdac_csda* v) {
    MDAC_SYM(MDAC_TRY {
        CSDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(common_env(l[i], ref(v)));
        assign(l[i], ref(v));
    } MDAC_CATCH)
}

mdac_status mdac_csdalist_push(mdac_csdalist* l_, const mdac_csda* v) {
    MDAC_SYM(MDAC_TRY {
        CSDAList& l = ref(l_);
        EnvGuard g(env_of(ref(v)));
        if (!l.empty()) common_env(l[0], ref(v));
        l.push_back(ref(v));
    } MDAC_CATCH)
}

// ---- cd_composition ----------------------------------------------------------------

mdac_status mdac_sdalist_compose_c(const mdac_sdalist* m, const mdac_csdalist* v, mdac_csdalist** out) {
    MDAC_SYM(return compose(ref(m), ref(v), out);)
}

mdac_status mdac_csdalist_compose(const mdac_csdalist* m, const mdac_csdalist* v, mdac_csdalist** out) {
    MDAC_SYM(return compose(ref(m), ref(v), out);)
}

mdac_status mdac_csdalist_compose_s(const mdac_csdalist* m, const mdac_sdalist* v, mdac_csdalist** out) {
    MDAC_SYM(return compose(ref(m), ref(v), out);)
}

}
