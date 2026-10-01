// capi_cnda.cpp — CNDA (std::complex<NDA>), CNDA lists and cd_composition in the
// C API (plan T4.1).
#include "common.h"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <string>

using namespace mdac;
using da::get_imag;
using da::get_real;
using Complex = std::complex<double>;

namespace {

da::DAEnv* env_of(const NDA& v) { return v.env_; }
da::DAEnv* env_of(const CNDA& v) { return get_real(v).env_; }

// The env shared by every DA argument.
template <class A, class... Rest>
da::DAEnv* common_env(const A& a, const Rest&... rest) {
    da::DAEnv* e = env_of(a);
    if (((env_of(rest) != e) || ...)) throw EnvError("DA vectors belong to different environments");
    return e;
}

// out = t part by part, keeping out's slots (a move assignment would swap them).
void assign(CNDA& out, const CNDA& t) {
    get_real(out) = get_real(t);
    get_imag(out) = get_imag(t);
}

// A new CNDA f() in the env of the DA arguments ds.
template <class F, class... D>
mdac_status make(mdac_cnda** out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ds...));
        *out = handle(new CNDA(f()));
    } MDAC_CATCH
}

// out = f() in out's own slots; out may be one of ds.
template <class F, class... D>
mdac_status into(mdac_cnda* out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ref(out), ds...));
        assign(ref(out), f());
    } MDAC_CATCH
}

template <class F, class... D>
mdac_status make_nda(mdac_nda** out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ds...));
        *out = handle(new NDA(f()));
    } MDAC_CATCH
}

// Every element of a list belongs to one env, so a list is freed in that env.
void check_joins(const CNDAList& l, const CNDA& v) {
    if (!l.empty()) common_env(l[0], v);
}

template <class L>
void same_env(da::DAEnv* e, const L& l) {
    for (const auto& v : l)
        if (env_of(v) != e) throw EnvError("DA vectors belong to different environments");
}

// cd_composition(m, v, out) as a new list; C++ only asserts the checks.
template <class M, class V>
mdac_status compose(const M& m, const V& v, mdac_cndalist** out) {
    MDAC_TRY {
        if (m.empty()) { *out = handle(new CNDAList()); return MDAC_OK; }
        da::DAEnv* e = env_of(m[0]);
        EnvGuard g(e);
        if (v.size() != e->layout().num_vars())
            throw std::invalid_argument("compose: the arguments must be nvars vectors");
        same_env(e, m);
        same_env(e, v);
        auto* o = new CNDAList(m.size());
        try {
            da::cd_composition(const_cast<M&>(m), const_cast<V&>(v), *o);
        } catch (...) {
            delete o;
            throw;
        }
        *out = handle(o);
    } MDAC_CATCH
}

} // namespace

extern "C" {

// ---- Lifecycle and inspection -------------------------------------------------

mdac_status mdac_cnda_new(const mdac_nda* re, const mdac_nda* im, mdac_cnda** out) {
    MDAC_TRY {
        EnvGuard g(im ? common_env(ref(re), ref(im)) : env_of(ref(re)));
        *out = handle(new CNDA(ref(re), im ? ref(im) : NDA(0.0)));
    } MDAC_CATCH
}

mdac_status mdac_cnda_new_z(mdac_env* e, double re, double im, mdac_cnda** out) {
    MDAC_TRY {
        EnvGuard g(env(e));
        *out = handle(new CNDA(NDA(re), NDA(im)));
    } MDAC_CATCH
}

mdac_status mdac_cnda_copy(const mdac_cnda* v, mdac_cnda** out) {
    return make(out, [&] { return ref(v); }, ref(v));
}

void mdac_cnda_free(mdac_cnda* v) { delete &ref(v); }

mdac_env* mdac_cnda_env(const mdac_cnda* v) { return handle(env_of(ref(v))); }

mdac_status mdac_cnda_real(const mdac_cnda* v, mdac_nda** out) {
    return make_nda(out, [&] { return get_real(ref(v)); }, ref(v));
}

mdac_status mdac_cnda_imag(const mdac_cnda* v, mdac_nda** out) {
    return make_nda(out, [&] { return get_imag(ref(v)); }, ref(v));
}

mdac_status mdac_cnda_set_real(mdac_cnda* v, const mdac_nda* x) {
    MDAC_TRY {
        EnvGuard g(common_env(ref(v), ref(x)));
        get_real(ref(v)) = ref(x);
    } MDAC_CATCH
}

mdac_status mdac_cnda_set_imag(mdac_cnda* v, const mdac_nda* x) {
    MDAC_TRY {
        EnvGuard g(common_env(ref(v), ref(x)));
        get_imag(ref(v)) = ref(x);
    } MDAC_CATCH
}

mdac_status mdac_cnda_to_string(const mdac_cnda* v, char* buf, size_t cap, size_t* n) {
    MDAC_TRY {
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
    } MDAC_CATCH
}

// ---- Arithmetic ------------------------------------------------------------------

#define MDAC_CNDA_PAIR(name, A, B, expr, ...)                                                \
    mdac_status name(A, B, mdac_cnda** out) { return make(out, [&] { return expr; }, __VA_ARGS__); } \
    mdac_status name##_into(mdac_cnda* out, A, B) {                                          \
        return into(out, [&] { return expr; }, __VA_ARGS__);                                 \
    }

#define MDAC_CNDA_ZPAIR(name, expr, ...)                                                     \
    mdac_status name(__VA_ARGS__, mdac_cnda** out) { return make(out, [&] { return expr; }, ref(a)); } \
    mdac_status name##_into(mdac_cnda* out, __VA_ARGS__) {                                   \
        return into(out, [&] { return expr; }, ref(a));                                      \
    }

#define MDAC_CNDA_BINOP(name, op)                                                            \
    MDAC_CNDA_PAIR(mdac_cnda_##name, const mdac_cnda* a, const mdac_cnda* b, ref(a) op ref(b), ref(a), ref(b)) \
    MDAC_CNDA_PAIR(mdac_cnda_##name##_n, const mdac_cnda* a, const mdac_nda* b, ref(a) op ref(b), ref(a), ref(b)) \
    MDAC_CNDA_PAIR(mdac_cnda_n##name, const mdac_nda* a, const mdac_cnda* b, ref(a) op ref(b), ref(a), ref(b)) \
    MDAC_CNDA_PAIR(mdac_cnda_##name##_d, const mdac_cnda* a, double x, ref(a) op x, ref(a)) \
    MDAC_CNDA_PAIR(mdac_cnda_d##name, double x, const mdac_cnda* a, x op ref(a), ref(a))   \
    MDAC_CNDA_ZPAIR(mdac_cnda_##name##_z, ref(a) op Complex(re, im), const mdac_cnda* a, double re, double im) \
    MDAC_CNDA_ZPAIR(mdac_cnda_z##name, Complex(re, im) op ref(a), double re, double im, const mdac_cnda* a) \
    MDAC_CNDA_ZPAIR(mdac_nda_##name##_z, ref(a) op Complex(re, im), const mdac_nda* a, double re, double im) \
    MDAC_CNDA_ZPAIR(mdac_nda_z##name, Complex(re, im) op ref(a), double re, double im, const mdac_nda* a)

MDAC_CNDA_BINOP(add, +)
MDAC_CNDA_BINOP(sub, -)
MDAC_CNDA_BINOP(mul, *)
MDAC_CNDA_BINOP(div, /)

MDAC_CNDA_PAIR(mdac_cnda_pow_i, const mdac_cnda* a, int n, da::pow(ref(a), n), ref(a))
MDAC_CNDA_PAIR(mdac_cnda_pow_d, const mdac_cnda* a, double x, da::pow(ref(a), x), ref(a))

mdac_status mdac_cnda_neg(const mdac_cnda* a, mdac_cnda** out) {
    return make(out, [&] { return -ref(a); }, ref(a));
}

mdac_status mdac_cnda_neg_into(mdac_cnda* out, const mdac_cnda* a) {
    return into(out, [&] { return -ref(a); }, ref(a));
}

// ---- Math functions --------------------------------------------------------------

#define MDAC_CNDA_FUNC(f)                                                                    \
    mdac_status mdac_cnda_##f(const mdac_cnda* a, mdac_cnda** out) {                        \
        return make(out, [&] { return da::f(ref(a)); }, ref(a));                             \
    }                                                                                        \
    mdac_status mdac_cnda_##f##_into(mdac_cnda* out, const mdac_cnda* a) {                  \
        return into(out, [&] { return da::f(ref(a)); }, ref(a));                             \
    }

MDAC_CNDA_FUNC(sqrt)
MDAC_CNDA_FUNC(exp)
MDAC_CNDA_FUNC(log)
MDAC_CNDA_FUNC(asin)
MDAC_CNDA_FUNC(acos)
MDAC_CNDA_FUNC(atan)
MDAC_CNDA_FUNC(asinh)
MDAC_CNDA_FUNC(acosh)
MDAC_CNDA_FUNC(atanh)

mdac_status mdac_cnda_abs(const mdac_cnda* a, double* out) {
    MDAC_TRY {
        EnvGuard g(env_of(ref(a)));
        *out = da::abs(ref(a));
    } MDAC_CATCH
}

// ---- Lists -------------------------------------------------------------------------

mdac_status mdac_cndalist_new(mdac_cndalist** out) {
    MDAC_TRY { *out = handle(new CNDAList()); } MDAC_CATCH
}

mdac_status mdac_cndalist_from(const mdac_cnda* const* vs, size_t n, mdac_cndalist** out) {
    MDAC_TRY {
        if (n == 0) { *out = handle(new CNDAList()); return MDAC_OK; }
        da::DAEnv* e = env_of(ref(vs[0]));
        EnvGuard g(e);
        for (size_t i = 1; i < n; ++i) common_env(ref(vs[0]), ref(vs[i]));
        auto* l = new CNDAList();
        try {
            l->reserve(n);
            for (size_t i = 0; i < n; ++i) l->push_back(ref(vs[i]));
        } catch (...) {
            delete l;
            throw;
        }
        *out = handle(l);
    } MDAC_CATCH
}

void mdac_cndalist_free(mdac_cndalist* l) { delete &ref(l); }

size_t mdac_cndalist_length(const mdac_cndalist* l) { return ref(l).size(); }

mdac_env* mdac_cndalist_env(const mdac_cndalist* l) {
    return ref(l).empty() ? nullptr : handle(env_of(ref(l)[0]));
}

mdac_status mdac_cndalist_get(const mdac_cndalist* l_, size_t i, mdac_cnda** out) {
    MDAC_TRY {
        const CNDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(env_of(l[i]));
        *out = handle(new CNDA(l[i]));
    } MDAC_CATCH
}

mdac_status mdac_cndalist_set(mdac_cndalist* l_, size_t i, const mdac_cnda* v) {
    MDAC_TRY {
        CNDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(common_env(l[i], ref(v)));
        assign(l[i], ref(v));
    } MDAC_CATCH
}

mdac_status mdac_cndalist_push(mdac_cndalist* l_, const mdac_cnda* v) {
    MDAC_TRY {
        CNDAList& l = ref(l_);
        EnvGuard g(env_of(ref(v)));
        check_joins(l, ref(v));
        l.push_back(ref(v));
    } MDAC_CATCH
}

// ---- cd_composition ----------------------------------------------------------------

mdac_status mdac_ndalist_compose_c(const mdac_ndalist* m, const mdac_cndalist* v, mdac_cndalist** out) {
    return compose(ref(m), ref(v), out);
}

mdac_status mdac_cndalist_compose(const mdac_cndalist* m, const mdac_cndalist* v, mdac_cndalist** out) {
    return compose(ref(m), ref(v), out);
}

mdac_status mdac_cndalist_compose_n(const mdac_cndalist* m, const mdac_ndalist* v, mdac_cndalist** out) {
    return compose(ref(m), ref(v), out);
}

}
