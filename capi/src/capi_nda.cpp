// capi_nda.cpp — NDA functions of the C API (plan T1.3-T1.5).
#include "common.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace mdac;

namespace {

// out = t, keeping out's slot (a move assignment would swap slots).
void assign(NDA& out, const NDA& t) { out = t; }

// The _into forms: out = a op b in out's own slot. When out is an operand,
// only out == a for + and - can work in place.
void check_into(const NDA& out, const NDA& a, const NDA& b) {
    check_same_env(out, a);
    check_same_env(out, b);
}

void add_into(NDA& out, const NDA& a, const NDA& b) {
    check_into(out, a, b);
    if (&out == &b && &out != &a) return assign(out, a + b);
    out = a;
    out += b;
}

void sub_into(NDA& out, const NDA& a, const NDA& b) {
    check_into(out, a, b);
    if (&out == &b && &out != &a) return assign(out, a - b);
    out = a;
    out -= b;
}

// The engine's three-slot kernels, which operator* and operator/ use too;
// NDA::operator*= would swap in a new slot.
void mul_into(NDA& out, const NDA& a, const NDA& b) {
    check_into(out, a, b);
    if (&out == &a || &out == &b) return assign(out, a * b);
    da::detail::ad_mult(out.env_->layout(), out.env_->pool<double>(), a.slot_, b.slot_, out.slot_);
}

void div_into(NDA& out, const NDA& a, const NDA& b) {
    check_into(out, a, b);
    if (&out == &a || &out == &b) return assign(out, a / b);
    da::detail::ad_div(out.env_->layout(), out.env_->pool<double>(), a.slot_, b.slot_, out.slot_);
}

void add_d_into(NDA& out, const NDA& a, double x) {
    check_same_env(out, a);
    out = a;
    out += x;
}

void sub_d_into(NDA& out, const NDA& a, double x) {
    check_same_env(out, a);
    out = a;
    out -= x;
}

void mul_d_into(NDA& out, const NDA& a, double x) {
    check_same_env(out, a);
    out = a;
    out *= x;
}

void div_d_into(NDA& out, const NDA& a, double x) {
    check_same_env(out, a);
    // The check of operator/(NDA, double).
    if (std::abs(x) < std::numeric_limits<double>::min())
        throw std::invalid_argument("da::operator/: divide by zero or a subnormal number");
    out = a;
    out /= x;
}

void dadd_into(NDA& out, double x, const NDA& a) { add_d_into(out, a, x); }

// x - a as operator-(double, NDA) computes it: -a + x.
void dsub_into(NDA& out, double x, const NDA& a) {
    check_same_env(out, a);
    out = a;
    out *= -1.0;
    out += x;
}

void dmul_into(NDA& out, double x, const NDA& a) { mul_d_into(out, a, x); }

void ddiv_into(NDA& out, double x, const NDA& a) {
    check_same_env(out, a);
    assign(out, x / a);
}

} // namespace

extern "C" {

// ---- Lifecycle and inspection (T1.3) --------------------------------------

mdac_status mdac_nda_new(mdac_env* e, double x, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(env(e));
        *out = handle(new NDA(x));
    } MDAC_CATCH
}

mdac_status mdac_nda_from_coeffs(mdac_env* e, const double* c, size_t n, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(env(e));
        if (n > env(e)->layout().full_len())
            throw std::invalid_argument("from_coeffs: more coefficients than full_length()");
        NDA* v = new NDA();
        da::Pool<double>& pool = v->env_->pool<double>();
        if (n > 0) std::memcpy(pool.slot(v->slot_), c, n * sizeof(double));
        pool.set_len(v->slot_, static_cast<unsigned>(std::max<size_t>(n, 1)));
        *out = handle(v);
    } MDAC_CATCH
}

mdac_status mdac_nda_var(mdac_env* e, unsigned i, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(env(e));
        if (i >= env(e)->layout().num_vars()) throw std::out_of_range("var: index out of range");
        *out = handle(new NDA(da::da_base(i)));
    } MDAC_CATCH
}

mdac_status mdac_nda_copy(const mdac_nda* v, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = handle(new NDA(ref(v)));
    } MDAC_CATCH
}

void mdac_nda_free(mdac_nda* v) { delete &ref(v); }

mdac_env* mdac_nda_env(const mdac_nda* v) { return handle(ref(v).env_); }

mdac_status mdac_nda_import(mdac_env* e, const mdac_nda* v, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = handle(new NDA(da::import_to(live(e), ref(v))));
    } MDAC_CATCH
}

#define MDAC_NDA_GET(name, T, expr)                                         \
    mdac_status mdac_nda_##name(const mdac_nda* v_, T* out) {               \
        MDAC_TRY { const NDA& v = ref(v_); EnvGuard g(v.env_); *out = (expr); } MDAC_CATCH \
    }

MDAC_NDA_GET(con, double, v.con())
MDAC_NDA_GET(length, size_t, v.length())
MDAC_NDA_GET(nterms, size_t, static_cast<size_t>(v.n_element()))
MDAC_NDA_GET(norm, double, v.norm())
MDAC_NDA_GET(abs, double, da::abs(v))

mdac_status mdac_nda_set_con(mdac_nda* v, double x) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).reset_const(x);
    } MDAC_CATCH
}

mdac_status mdac_nda_coeffs(const mdac_nda* v_, double* buf, size_t cap, size_t* n) {
    MDAC_TRY {
        const NDA& v = ref(v_);
        EnvGuard g(v.env_);
        *n = v.length();
        std::memcpy(buf, v.env_->pool<double>().slot(v.slot_), std::min(cap, *n) * sizeof(double));
    } MDAC_CATCH
}

mdac_status mdac_nda_coeff(const mdac_nda* v, const int* exps, size_t k, double* out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = const_cast<NDA&>(ref(v)).element(exponents(exps, k));
    } MDAC_CATCH
}

mdac_status mdac_nda_set_coeff(mdac_nda* v, const int* exps, size_t k, double x) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).set_element(exponents(exps, k), x);
    } MDAC_CATCH
}

mdac_status mdac_nda_index_term(const mdac_nda* v_, size_t i, int* exps, double* out) {
    MDAC_TRY {
        const NDA& v = ref(v_);
        EnvGuard g(v.env_);
        if (i >= v.env_->layout().full_len()) throw std::out_of_range("index_term: index out of range");
        std::vector<unsigned> c;
        v.element(static_cast<unsigned>(i), c, *out);
        std::copy(c.begin(), c.end(), exps);
    } MDAC_CATCH
}

mdac_status mdac_nda_iszero(const mdac_nda* v, double eps, int* out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = ref(v).iszero(eps) ? 1 : 0;
    } MDAC_CATCH
}

mdac_status mdac_nda_clean(mdac_nda* v, double eps) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).clean(eps);
    } MDAC_CATCH
}

mdac_status mdac_nda_reset(mdac_nda* v) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).reset();
    } MDAC_CATCH
}

mdac_status mdac_nda_to_string(const mdac_nda* v, char* buf, size_t cap, size_t* n) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
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

// ---- Arithmetic (T1.4) ------------------------------------------------------

#define MDAC_NDA_BINOP(name, op)                                                            \
    mdac_status mdac_nda_##name(const mdac_nda* a, const mdac_nda* b, mdac_nda** out) {     \
        MDAC_TRY { EnvGuard g(ref(a).env_); *out = handle(new NDA(ref(a) op ref(b))); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_##name##_d(const mdac_nda* a, double x, mdac_nda** out) {          \
        MDAC_TRY { EnvGuard g(ref(a).env_); *out = handle(new NDA(ref(a) op x)); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_d##name(double x, const mdac_nda* a, mdac_nda** out) {             \
        MDAC_TRY { EnvGuard g(ref(a).env_); *out = handle(new NDA(x op ref(a))); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_##name##_into(mdac_nda* out, const mdac_nda* a, const mdac_nda* b) { \
        MDAC_TRY { EnvGuard g(ref(a).env_); name##_into(ref(out), ref(a), ref(b)); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_##name##_d_into(mdac_nda* out, const mdac_nda* a, double x) {      \
        MDAC_TRY { EnvGuard g(ref(a).env_); name##_d_into(ref(out), ref(a), x); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_d##name##_into(mdac_nda* out, double x, const mdac_nda* a) {       \
        MDAC_TRY { EnvGuard g(ref(a).env_); d##name##_into(ref(out), x, ref(a)); } MDAC_CATCH \
    }

MDAC_NDA_BINOP(add, +)
MDAC_NDA_BINOP(sub, -)
MDAC_NDA_BINOP(mul, *)
MDAC_NDA_BINOP(div, /)

// f(a, args...) in both forms; the _into form assigns the result to out.
#define MDAC_NDA_UNARY(name, expr, ...)                                                     \
    mdac_status mdac_nda_##name(const mdac_nda* a_ __VA_ARGS__, mdac_nda** out) {           \
        MDAC_TRY { const NDA& a = ref(a_); EnvGuard g(a.env_); *out = handle(new NDA(expr)); } MDAC_CATCH \
    }                                                                                       \
    mdac_status mdac_nda_##name##_into(mdac_nda* out, const mdac_nda* a_ __VA_ARGS__) {     \
        MDAC_TRY {                                                                          \
            const NDA& a = ref(a_);                                                         \
            EnvGuard g(a.env_);                                                             \
            check_same_env(ref(out), a);                                                          \
            assign(ref(out), expr);                                                         \
        } MDAC_CATCH                                                                        \
    }

MDAC_NDA_UNARY(neg, -a)
MDAC_NDA_UNARY(pow_i, da::pow(a, n), , int n)
MDAC_NDA_UNARY(pow_d, da::pow(a, x), , double x)

// ---- Math functions (T1.5) --------------------------------------------------

#define MDAC_NDA_FUNC(f) MDAC_NDA_UNARY(f, da::f(a))

MDAC_NDA_FUNC(sqrt)
MDAC_NDA_FUNC(exp)
MDAC_NDA_FUNC(log)
MDAC_NDA_FUNC(sin)
MDAC_NDA_FUNC(cos)
MDAC_NDA_FUNC(tan)
MDAC_NDA_FUNC(asin)
MDAC_NDA_FUNC(acos)
MDAC_NDA_FUNC(atan)
MDAC_NDA_FUNC(sinh)
MDAC_NDA_FUNC(cosh)
MDAC_NDA_FUNC(tanh)
MDAC_NDA_FUNC(asinh)
MDAC_NDA_FUNC(acosh)
MDAC_NDA_FUNC(atanh)
MDAC_NDA_FUNC(erf)

}
