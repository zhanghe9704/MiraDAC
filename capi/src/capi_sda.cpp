// capi_sda.cpp — Expr (SymEngine::Expression), SDA, SDA lists and the SDA
// algorithms of the C API (plan T5.1, T5.2). Without DA_WITH_SYMBOLIC every
// function is a stub (MDAC_SYM, MDAC_SYM_ELSE).
//
// NDA functions with no SDA form: from_coeffs, norm, abs and eval at a point
// (the symbolic kernels count non-zero coefficients as 0 or need numbers),
// the eps argument of iszero and clean (the symbolic kernels test with
// is_zero()), asinh, acosh, atanh, inv_map and evaluate_map (no SDA version
// in C++).
#include "common.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef DA_WITH_SYMBOLIC
#include <symengine/parser.h>
#include <symengine/visitor.h>

using namespace mdac;

namespace {

mdac_expr* new_expr(Expression x) { return handle(new Expression(std::move(x))); }

// The algorithms take non-const lists but do not change their inputs.
SDAList& mut(const SDAList& l) { return const_cast<SDAList&>(l); }

// The operator<< text of v (size-query convention, as mdac_nda_to_string).
template <class T>
void to_string(const T& v, char* buf, size_t cap, size_t* n) {
    std::ostringstream os;
    os << v;
    const std::string s = os.str();
    *n = s.size() + 1;
    if (cap > 0) {
        const size_t m = std::min(cap - 1, s.size());
        std::memcpy(buf, s.data(), m);
        buf[m] = '\0';
    }
}

SymEngine::map_basic_basic to_map(size_t n, const mdac_expr* const* keys,
                                  const mdac_expr* const* vals) {
    SymEngine::map_basic_basic m;
    for (size_t i = 0; i < n; ++i) m[ref(keys[i]).get_basic()] = ref(vals[i]).get_basic();
    return m;
}

da::DAEnv* env_of(const NDA& v) { return v.env_; }
da::DAEnv* env_of(const SDA& v) { return v.env_; }

// The env shared by every DA argument.
template <class A, class... Rest>
da::DAEnv* common_env(const A& a, const Rest&... rest) {
    da::DAEnv* e = env_of(a);
    if (((env_of(rest) != e) || ...)) throw EnvError("DA vectors belong to different environments");
    return e;
}

template <class F>
mdac_status make_expr(mdac_expr** out, F f) {
    MDAC_TRY { *out = new_expr(f()); } MDAC_CATCH
}

// A new SDA f() in the env of the DA arguments ds.
template <class F, class... D>
mdac_status make(mdac_sda** out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ds...));
        *out = handle(new SDA(f()));
    } MDAC_CATCH
}

// out = f() in out's own slot (a move assignment would swap slots); out may
// be one of ds.
template <class F, class... D>
mdac_status into(mdac_sda* out, F f, const D&... ds) {
    MDAC_TRY {
        EnvGuard g(common_env(ref(out), ds...));
        const SDA t = f();
        ref(out) = t;
    } MDAC_CATCH
}

// A new SDA holding f(coefficient) for every coefficient of v.
template <class F>
mdac_status map_coeffs(const SDA& v, mdac_sda** out, F f) {
    MDAC_TRY {
        EnvGuard g(v.env_);
        auto r = std::make_unique<SDA>(v);
        da::Pool<Expression>& pool = r->env_->pool<Expression>();
        Expression* p = pool.slot(r->slot_);
        for (unsigned i = 0, n = pool.len(r->slot_); i < n; ++i) f(p[i]);
        *out = handle(r.release());
    } MDAC_CATCH
}

void check_base(const SDA& v, unsigned id) {
    if (id >= v.env_->layout().num_vars()) throw std::out_of_range("base id out of range");
}

void same_env(const SDA& first, const SDAList& l) {
    for (const SDA& v : l) check_same_env(first, v);
}

std::vector<unsigned> check_ids(const SDA& first, const unsigned* ids, size_t k, const SDAList& v) {
    if (k != v.size()) throw std::invalid_argument("base ids and vectors must have the same length");
    std::vector<unsigned> b(ids, ids + k);
    if (std::set<unsigned>(b.begin(), b.end()).size() != b.size())
        throw std::invalid_argument("duplicate base id");
    for (unsigned id : b) check_base(first, id);
    same_env(first, v);
    return b;
}

// A new list of n zero vectors in the current env.
std::unique_ptr<SDAList> zeros(size_t n) {
    auto l = std::make_unique<SDAList>();
    l->reserve(n);
    for (size_t i = 0; i < n; ++i) l->emplace_back();
    return l;
}

// A new SDA computed into o by f(o).
template <class F>
mdac_status make_into(const SDA& iv, mdac_sda** out, F f) {
    MDAC_TRY {
        EnvGuard g(iv.env_);
        auto o = std::make_unique<SDA>();
        f(*o);
        *out = handle(o.release());
    } MDAC_CATCH
}

} // namespace
#endif // DA_WITH_SYMBOLIC

extern "C" {

// ---- Expr (T5.1) ------------------------------------------------------------

mdac_status mdac_expr_new_d(double x, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return Expression(x); });)
}

mdac_status mdac_expr_new_i(int64_t x, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return Expression(SymEngine::integer(x)); });)
}

mdac_status mdac_expr_parse(const char* s, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return Expression(SymEngine::parse(s)); });)
}

mdac_status mdac_expr_symbol(const char* name, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return Expression(SymEngine::symbol(name)); });)
}

mdac_status mdac_expr_copy(const mdac_expr* x, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return ref(x); });)
}

void mdac_expr_free(mdac_expr* x) { MDAC_SYM_ELSE((void)x;, delete &ref(x);) }

mdac_status mdac_expr_to_string(const mdac_expr* x, char* buf, size_t cap, size_t* n) {
    MDAC_SYM(MDAC_TRY { to_string(ref(x), buf, cap, n); } MDAC_CATCH)
}

mdac_status mdac_expr_to_double(const mdac_expr* x, double* out) {
    MDAC_SYM(MDAC_TRY {
        if (!SymEngine::free_symbols(*ref(x).get_basic()).empty())
            throw std::invalid_argument("cannot convert an expression with free symbols to a number");
        *out = SymEngine::eval_double(*ref(x).get_basic());
    } MDAC_CATCH)
}

#define MDAC_EXPR_BINOP(name, expr)                                                          \
    mdac_status mdac_expr_##name(const mdac_expr* a_, const mdac_expr* b_, mdac_expr** out) { \
        MDAC_SYM(return make_expr(out, [&] { const Expression &a = ref(a_), &b = ref(b_); return expr; });) \
    }                                                                                        \
    mdac_status mdac_expr_##name##_d(const mdac_expr* a_, double x, mdac_expr** out) {       \
        MDAC_SYM(return make_expr(out, [&] { const Expression &a = ref(a_), b(x); return expr; });) \
    }                                                                                        \
    mdac_status mdac_expr_d##name(double x, const mdac_expr* b_, mdac_expr** out) {          \
        MDAC_SYM(return make_expr(out, [&] { const Expression a(x), &b = ref(b_); return expr; });) \
    }

MDAC_EXPR_BINOP(add, a + b)
MDAC_EXPR_BINOP(sub, a - b)
MDAC_EXPR_BINOP(mul, a * b)
MDAC_EXPR_BINOP(div, a / b)
MDAC_EXPR_BINOP(pow, SymEngine::pow(a, b))

mdac_status mdac_expr_neg(const mdac_expr* a, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return -ref(a); });)
}

mdac_status mdac_expr_eq(const mdac_expr* a, const mdac_expr* b, int* out) {
    MDAC_SYM(MDAC_TRY { *out = ref(a) == ref(b); } MDAC_CATCH)
}

mdac_status mdac_expr_hash(const mdac_expr* a, uint64_t* out) {
    MDAC_SYM(MDAC_TRY { *out = static_cast<uint64_t>(ref(a).get_basic()->hash()); } MDAC_CATCH)
}

mdac_status mdac_expr_is_zero(const mdac_expr* x, int* out) {
    MDAC_SYM(MDAC_TRY { *out = da::is_zero(ref(x)); } MDAC_CATCH)
}

mdac_status mdac_expr_subs(const mdac_expr* x, size_t n, const mdac_expr* const* keys,
                           const mdac_expr* const* vals, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return ref(x).subs(to_map(n, keys, vals)); });)
}

mdac_status mdac_expr_expand(const mdac_expr* x, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] { return Expression(SymEngine::expand(ref(x))); });)
}

mdac_status mdac_expr_diff(const mdac_expr* x, const mdac_expr* sym, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] {
        const auto& s = ref(sym).get_basic();
        if (!SymEngine::is_a<SymEngine::Symbol>(*s))
            throw std::invalid_argument("diff: the variable must be a symbol");
        return ref(x).diff(SymEngine::rcp_static_cast<const SymEngine::Symbol>(s));
    });)
}

mdac_status mdac_expr_simplify(const mdac_expr* x, mdac_expr** out) {
    MDAC_SYM(return make_expr(out, [&] {
        Expression r = ref(x);
        da::simplified_expr(r);
        return r;
    });)
}

mdac_status mdac_expr_free_symbols(const mdac_expr* x, mdac_expr** buf, size_t cap, size_t* n) {
    MDAC_SYM(MDAC_TRY {
        const SymEngine::set_basic s = SymEngine::free_symbols(*ref(x).get_basic());
        std::vector<Expression> syms(s.begin(), s.end());
        const size_t m = std::min(cap, syms.size());
        std::vector<std::unique_ptr<Expression>> made;
        for (size_t i = 0; i < m; ++i) made.push_back(std::make_unique<Expression>(syms[i]));
        for (size_t i = 0; i < m; ++i) buf[i] = handle(made[i].release());
        *n = syms.size();
    } MDAC_CATCH)
}

// ---- SDA: lifecycle and inspection (T5.2) ----------------------------------

mdac_status mdac_sda_new(mdac_env* e, const mdac_expr* x, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env(e));
        *out = handle(x ? new SDA(ref(x)) : new SDA());
    } MDAC_CATCH)
}

mdac_status mdac_sda_new_d(mdac_env* e, double x, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env(e));
        *out = handle(new SDA(x));
    } MDAC_CATCH)
}

mdac_status mdac_sda_var(mdac_env* e, unsigned i, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(env(e));
        if (i >= env(e)->layout().num_vars()) throw std::out_of_range("var: index out of range");
        *out = handle(new SDA(da::promote(da::da_base(i))));
    } MDAC_CATCH)
}

mdac_status mdac_sda_promote(const mdac_nda* v, mdac_sda** out) {
    MDAC_SYM(return make(out, [&] { return da::promote(ref(v)); }, ref(v));)
}

mdac_status mdac_sda_copy(const mdac_sda* v, mdac_sda** out) {
    MDAC_SYM(return make(out, [&] { return SDA(ref(v)); }, ref(v));)
}

void mdac_sda_free(mdac_sda* v) { MDAC_SYM_ELSE((void)v;, delete &ref(v);) }

mdac_env* mdac_sda_env(const mdac_sda* v) {
    MDAC_SYM_ELSE((void)v; return nullptr;, return mdac::handle(ref(v).env_);)
}

mdac_status mdac_sda_import(mdac_env* e, const mdac_sda* v, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = handle(new SDA(da::import_to(live(e), ref(v))));
    } MDAC_CATCH)
}

mdac_status mdac_nda_promote_to(mdac_env* e, const mdac_nda* v, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = handle(new SDA(da::promote_to(live(e), ref(v))));
    } MDAC_CATCH)
}

#define MDAC_SDA_GET(name, T, expr)                                           \
    mdac_status mdac_sda_##name(const mdac_sda* v_, T* out) {                 \
        MDAC_SYM(MDAC_TRY { const SDA& v = ref(v_); EnvGuard g(v.env_); *out = (expr); } MDAC_CATCH) \
    }

MDAC_SDA_GET(con, mdac_expr*, new_expr(v.con()))
MDAC_SDA_GET(length, size_t, v.length())
MDAC_SDA_GET(nterms, size_t, static_cast<size_t>(v.n_element()))
MDAC_SDA_GET(iszero, int, v.iszero() ? 1 : 0)

mdac_status mdac_sda_set_con(mdac_sda* v, const mdac_expr* x) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).reset_const_expr(ref(x));
    } MDAC_CATCH)
}

mdac_status mdac_sda_coeffs(const mdac_sda* v_, mdac_expr** buf, size_t cap, size_t* n) {
    MDAC_SYM(MDAC_TRY {
        const SDA& v = ref(v_);
        EnvGuard g(v.env_);
        const Expression* p = v.env_->pool<Expression>().slot(v.slot_);
        const size_t len = v.length(), m = std::min(cap, len);
        std::vector<std::unique_ptr<Expression>> made;
        for (size_t i = 0; i < m; ++i) made.push_back(std::make_unique<Expression>(p[i]));
        for (size_t i = 0; i < m; ++i) buf[i] = handle(made[i].release());
        *n = len;
    } MDAC_CATCH)
}

mdac_status mdac_sda_coeff(const mdac_sda* v, const int* exps, size_t k, mdac_expr** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        *out = new_expr(const_cast<SDA&>(ref(v)).element_expr(exponents(exps, k)));
    } MDAC_CATCH)
}

mdac_status mdac_sda_set_coeff(mdac_sda* v, const int* exps, size_t k, const mdac_expr* x) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        ref(v).set_element(exponents(exps, k), ref(x));
    } MDAC_CATCH)
}

mdac_status mdac_sda_index_term(const mdac_sda* v_, size_t i, int* exps, mdac_expr** out) {
    MDAC_SYM(MDAC_TRY {
        const SDA& v = ref(v_);
        EnvGuard g(v.env_);
        da::Layout& l = v.env_->layout();
        if (i >= l.full_len()) throw std::out_of_range("index_term: index out of range");
        std::vector<unsigned> c(l.num_vars());
        Expression value;
        da::detail::ad_elem(l, v.env_->pool<Expression>(), v.slot_, static_cast<unsigned>(i) + 1,
                            c.data(), value);
        std::copy(c.begin(), c.end(), exps);
        *out = new_expr(value);
    } MDAC_CATCH)
}

mdac_status mdac_sda_clean(mdac_sda* v) {
    MDAC_SYM(MDAC_TRY { EnvGuard g(ref(v).env_); ref(v).clean(); } MDAC_CATCH)
}

mdac_status mdac_sda_reset(mdac_sda* v) {
    MDAC_SYM(MDAC_TRY { EnvGuard g(ref(v).env_); ref(v).reset(); } MDAC_CATCH)
}

mdac_status mdac_sda_to_string(const mdac_sda* v, char* buf, size_t cap, size_t* n) {
    MDAC_SYM(MDAC_TRY { EnvGuard g(ref(v).env_); to_string(ref(v), buf, cap, n); } MDAC_CATCH)
}

// ---- SDA: arithmetic --------------------------------------------------------
// One allocating and one _into form per operand pair (MDAC_SDA_OP2): A and B
// are the operands as C++ objects, ds the DA ones (which give the env).

#define MDAC_UNPAREN(...) __VA_ARGS__

#define MDAC_SDA_BINOP(op, op_)                                                              \
    MDAC_SDA_OP2(mdac_sda_##op, ref(a), ref(b), (ref(a), ref(b)), op_, const mdac_sda* a, const mdac_sda* b) \
    MDAC_SDA_OP2(mdac_sda_##op##_n, ref(a), ref(b), (ref(a), ref(b)), op_, const mdac_sda* a, const mdac_nda* b) \
    MDAC_SDA_OP2(mdac_sda_n##op, ref(a), ref(b), (ref(a), ref(b)), op_, const mdac_nda* a, const mdac_sda* b) \
    MDAC_SDA_OP2(mdac_sda_##op##_e, ref(a), ref(x), (ref(a)), op_, const mdac_sda* a, const mdac_expr* x) \
    MDAC_SDA_OP2(mdac_sda_e##op, ref(x), ref(a), (ref(a)), op_, const mdac_expr* x, const mdac_sda* a) \
    MDAC_SDA_OP2(mdac_sda_##op##_d, ref(a), x, (ref(a)), op_, const mdac_sda* a, double x)  \
    MDAC_SDA_OP2(mdac_sda_d##op, x, ref(a), (ref(a)), op_, double x, const mdac_sda* a)     \
    MDAC_SDA_OP2(mdac_nda_##op##_e, ref(a), ref(x), (ref(a)), op_, const mdac_nda* a, const mdac_expr* x) \
    MDAC_SDA_OP2(mdac_nda_e##op, ref(x), ref(a), (ref(a)), op_, const mdac_expr* x, const mdac_nda* a)

#define MDAC_SDA_OP2(name, A, B, ds, op_, ...)                                               \
    mdac_status name(__VA_ARGS__, mdac_sda** out) {                                          \
        MDAC_SYM(return make(out, [&] { return A op_ B; }, MDAC_UNPAREN ds);)                \
    }                                                                                        \
    mdac_status name##_into(mdac_sda* out, __VA_ARGS__) {                                    \
        MDAC_SYM(return into(out, [&] { return A op_ B; }, MDAC_UNPAREN ds);)                \
    }

MDAC_SDA_BINOP(add, +)
MDAC_SDA_BINOP(sub, -)
MDAC_SDA_BINOP(mul, *)
MDAC_SDA_BINOP(div, /)

// f(a, args...) in both forms.
#define MDAC_SDA_UNARY(name, expr, ...)                                                      \
    mdac_status mdac_sda_##name(const mdac_sda* a_ __VA_ARGS__, mdac_sda** out) {            \
        MDAC_SYM(const SDA& a = ref(a_); return make(out, [&] { return expr; }, a);)         \
    }                                                                                        \
    mdac_status mdac_sda_##name##_into(mdac_sda* out, const mdac_sda* a_ __VA_ARGS__) {      \
        MDAC_SYM(const SDA& a = ref(a_); return into(out, [&] { return expr; }, a);)         \
    }

MDAC_SDA_UNARY(neg, -a)
MDAC_SDA_UNARY(pow_i, da::pow(a, n), , int n)
MDAC_SDA_UNARY(pow_d, da::pow(a, x), , double x)

#define MDAC_SDA_FUNC(f) MDAC_SDA_UNARY(f, da::f(a))

MDAC_SDA_FUNC(sqrt)
MDAC_SDA_FUNC(exp)
MDAC_SDA_FUNC(log)
MDAC_SDA_FUNC(sin)
MDAC_SDA_FUNC(cos)
MDAC_SDA_FUNC(tan)
MDAC_SDA_FUNC(asin)
MDAC_SDA_FUNC(acos)
MDAC_SDA_FUNC(atan)
MDAC_SDA_FUNC(sinh)
MDAC_SDA_FUNC(cosh)
MDAC_SDA_FUNC(tanh)
MDAC_SDA_FUNC(erf)

// ---- SDA: per-coefficient operations and evaluation -------------------------

mdac_status mdac_sda_simplify(const mdac_sda* v, mdac_sda** out) {
    MDAC_SYM(return map_coeffs(ref(v), out, [](Expression& e) {
        if (!da::is_zero(e)) da::simplified_expr(e);
    });)
}

mdac_status mdac_sda_expand(const mdac_sda* v, mdac_sda** out) {
    MDAC_SYM(return map_coeffs(ref(v), out, [](Expression& e) { e = SymEngine::expand(e); });)
}

mdac_status mdac_sda_subs(const mdac_sda* v, size_t n, const mdac_expr* const* keys,
                          const mdac_expr* const* vals, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        const SymEngine::map_basic_basic m = to_map(n, keys, vals);
        return map_coeffs(ref(v), out, [&](Expression& e) { e = e.subs(m); });
    } MDAC_CATCH)
}

// The vector overload of da::evaluate (interop.h), as the Python binding uses.
mdac_status mdac_sda_evaluate(const mdac_sda* v, size_t n, const mdac_expr* const* keys,
                              const double* vals, mdac_nda** out) {
    MDAC_SYM(MDAC_TRY {
        EnvGuard g(ref(v).env_);
        std::vector<SymEngine::RCP<const SymEngine::Basic>> syms;
        for (size_t i = 0; i < n; ++i) syms.push_back(ref(keys[i]).get_basic());
        const std::vector<double> x(vals, vals + n);
        *out = mdac::handle(new NDA(da::evaluate(ref(v), syms, x)));
    } MDAC_CATCH)
}

// ---- SDA lists ----------------------------------------------------------------

mdac_status mdac_sdalist_new(mdac_sdalist** out) {
    MDAC_SYM(MDAC_TRY { *out = handle(new SDAList()); } MDAC_CATCH)
}

mdac_status mdac_sdalist_from(const mdac_sda* const* vs, size_t n, mdac_sdalist** out) {
    MDAC_SYM(MDAC_TRY {
        auto l = std::make_unique<SDAList>();
        if (n > 0) {
            EnvGuard g(ref(vs[0]).env_);
            for (size_t i = 1; i < n; ++i) check_same_env(ref(vs[0]), ref(vs[i]));
            l->reserve(n);
            for (size_t i = 0; i < n; ++i) l->push_back(ref(vs[i]));
        }
        *out = handle(l.release());
    } MDAC_CATCH)
}

void mdac_sdalist_free(mdac_sdalist* l) { MDAC_SYM_ELSE((void)l;, delete &ref(l);) }

size_t mdac_sdalist_length(const mdac_sdalist* l) {
    MDAC_SYM_ELSE((void)l; return 0;, return ref(l).size();)
}

mdac_env* mdac_sdalist_env(const mdac_sdalist* l) {
    MDAC_SYM_ELSE((void)l; return nullptr;,
                  return ref(l).empty() ? nullptr : mdac::handle(ref(l)[0].env_);)
}

mdac_status mdac_sdalist_get(const mdac_sdalist* l_, size_t i, mdac_sda** out) {
    MDAC_SYM(MDAC_TRY {
        const SDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(l[i].env_);
        *out = handle(new SDA(l[i]));
    } MDAC_CATCH)
}

mdac_status mdac_sdalist_set(mdac_sdalist* l_, size_t i, const mdac_sda* v) {
    MDAC_SYM(MDAC_TRY {
        SDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(ref(v).env_);
        check_same_env(l[i], ref(v));
        l[i] = ref(v);
    } MDAC_CATCH)
}

mdac_status mdac_sdalist_push(mdac_sdalist* l_, const mdac_sda* v) {
    MDAC_SYM(MDAC_TRY {
        SDAList& l = ref(l_);
        EnvGuard g(ref(v).env_);
        if (!l.empty()) check_same_env(l[0], ref(v));
        l.push_back(ref(v));
    } MDAC_CATCH)
}

// ---- SDA algorithms -----------------------------------------------------------

mdac_status mdac_sda_der(const mdac_sda* v, unsigned i, mdac_sda** out) {
    MDAC_SYM(return make(out, [&] { check_base(ref(v), i); return da::da_der(ref(v), i); }, ref(v));)
}

mdac_status mdac_sda_integ(const mdac_sda* v, unsigned i, mdac_sda** out) {
    MDAC_SYM(return make(out, [&] { check_base(ref(v), i); return da::da_int(ref(v), i); }, ref(v));)
}

mdac_status mdac_sda_substitute_d(const mdac_sda* iv, unsigned i, double x, mdac_sda** out) {
    MDAC_SYM(return make_into(ref(iv), out, [&](SDA& o) {
        check_base(ref(iv), i);
        da::da_substitute_const(ref(iv), i, x, o);
    });)
}

mdac_status mdac_sda_substitute(const mdac_sda* iv, unsigned i, const mdac_sda* v, mdac_sda** out) {
    MDAC_SYM(return make_into(ref(iv), out, [&](SDA& o) {
        check_same_env(ref(iv), ref(v));
        check_base(ref(iv), i);
        da::da_substitute(ref(iv), i, ref(v), o);
    });)
}

mdac_status mdac_sda_substitute_multi(const mdac_sda* iv, const unsigned* ids, size_t k,
                                      const mdac_sdalist* v, mdac_sda** out) {
    MDAC_SYM(return make_into(ref(iv), out, [&](SDA& o) {
        std::vector<unsigned> b = check_ids(ref(iv), ids, k, ref(v));
        da::da_substitute(ref(iv), b, mut(ref(v)), o);
    });)
}

mdac_status mdac_sdalist_substitute(const mdac_sdalist* m_, const unsigned* ids, size_t k,
                                    const mdac_sdalist* v, mdac_sdalist** out) {
    MDAC_SYM(MDAC_TRY {
        const SDAList& m = ref(m_);
        if (m.empty()) { *out = handle(new SDAList()); return MDAC_OK; }
        EnvGuard g(m[0].env_);
        std::vector<unsigned> b = check_ids(m[0], ids, k, ref(v));
        same_env(m[0], m);
        auto o = zeros(m.size());
        da::da_substitute(mut(m), b, mut(ref(v)), *o);
        *out = handle(o.release());
    } MDAC_CATCH)
}

mdac_status mdac_sdalist_compose(const mdac_sdalist* m_, const mdac_sdalist* v_, mdac_sdalist** out) {
    MDAC_SYM(MDAC_TRY {
        const SDAList &m = ref(m_), &v = ref(v_);
        if (m.empty()) { *out = handle(new SDAList()); return MDAC_OK; }
        EnvGuard g(m[0].env_);
        if (v.size() != m[0].env_->layout().num_vars())
            throw std::invalid_argument("compose: the arguments must be nvars vectors");
        same_env(m[0], m);
        same_env(m[0], v);
        auto o = zeros(m.size());
        da::da_composition(mut(m), mut(v), *o);
        *out = handle(o.release());
    } MDAC_CATCH)
}

mdac_status mdac_sdalist_compose_d(const mdac_sdalist* m_, const double* pt, size_t n,
                                   mdac_expr** out) {
    MDAC_SYM(MDAC_TRY {
        const SDAList& m = ref(m_);
        if (m.empty()) return MDAC_OK;
        EnvGuard g(m[0].env_);
        if (n != m[0].env_->layout().num_vars())
            throw std::invalid_argument("the point must have nvars coordinates");
        same_env(m[0], m);
        std::vector<double> p(pt, pt + n);
        std::vector<Expression> o;
        da::da_composition(mut(m), p, o);
        std::vector<std::unique_ptr<Expression>> made;
        for (const Expression& e : o) made.push_back(std::make_unique<Expression>(e));
        for (size_t i = 0; i < made.size(); ++i) out[i] = handle(made[i].release());
    } MDAC_CATCH)
}

}
