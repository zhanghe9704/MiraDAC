// bind_sda.cpp — SDA (DAVector<SymEngine::Expression>), SDAList, promote(),
// svar(), evaluate(), the SDA math functions and algorithms (plan T5.1-T5.6).
//
// NDA methods not bound for SDA: norm() and weighted_norm() (the kernels
// count every non-zero symbolic coefficient as 0), from_coeffs() (use
// promote(NDA.from_coeffs(a))), __call__ (use evaluate() and call the NDA),
// and the eps argument of iszero() and clean() (the symbolic kernels ignore
// it and test with is_zero()). No SDA versions exist in C++ of asinh, acosh,
// atanh, abs, inv_map or the file helpers.
#include "arith.h"
#include "se_bridge.h"

#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <symengine/basic.h>

#include <set>
#include <sstream>

using SymEngine::Expression;
using se_bridge::from_int;
using se_bridge::to_expr;

namespace {

constexpr auto symengine_operand = [](const se_bridge::SymengineBasic& h) {
    return se_bridge::from_symengine(h);
};
constexpr auto int_operand = [](const nb::int_& i) { return from_int(i); };

// Binds SDA (op) S and S (op) SDA, and SDA (op)= S, converting S to an
// Expression with conv first.
template <class S, class Conv>
void exprops(nb::class_<SDA>& c, Conv conv) {
    using namespace detail_arith;
    auto both = [&](const char* name, const char* rname, auto op) {
        c.def(name, [conv, op](const SDA& a, const S& x) {
            Expression e = conv(x);
            EnvGuard g(a.env_);
            return make_da<SDA>([&] { return op(a, e); });
        }, nb::is_operator(), op_doc(name));
        c.def(rname, [conv, op](const SDA& a, const S& x) {
            Expression e = conv(x);
            EnvGuard g(a.env_);
            return make_da<SDA>([&] { return op(e, a); });
        }, nb::is_operator(), op_doc(rname));
    };
    both("__add__", "__radd__", add);
    both("__sub__", "__rsub__", sub);
    both("__mul__", "__rmul__", mul);
    both("__truediv__", "__rtruediv__", divide);
    inplace<S>(c, "__iadd__", [conv](SDA& a, const S& x) { a += conv(x); });
    inplace<S>(c, "__isub__", [conv](SDA& a, const S& x) { a -= conv(x); });
    inplace<S>(c, "__imul__", [conv](SDA& a, const S& x) { a *= conv(x); });
    inplace<S>(c, "__itruediv__", [conv](SDA& a, const S& x) { a /= conv(x); });
}

template <SDA (*F)(const SDA&)>
auto unary(const SDA& v) {
    EnvGuard g(v.env_);
    return make_da<SDA>([&] { return F(v); });
}

// A new SDA holding f(coefficient) for every coefficient of v (T5.6).
template <class F>
auto map_coeffs(const SDA& v, F f) {
    EnvGuard g(v.env_);
    auto o = make_da<SDA>([&] { return SDA(v); });
    SDA& r = nb::cast<SDA&>(o);
    da::Pool<Expression>& pool = r.env_->pool<Expression>();
    Expression* p = pool.slot(r.slot_);
    for (unsigned i = 0, n = pool.len(r.slot_); i < n; ++i) f(p[i]);
    return o;
}

SymEngine::map_basic_basic to_map(const nb::typed<nb::dict, nb::any, nb::any>& d) {
    SymEngine::map_basic_basic map;
    for (auto [k, v] : d) map[to_expr(k).get_basic()] = to_expr(v).get_basic();
    return map;
}

void check_multi(const SDA& first, const std::vector<unsigned>& ids, const SDAList& v) {
    if (ids.size() != v.size())
        throw std::invalid_argument("base_id and v must have the same length");
    if (std::set<unsigned>(ids.begin(), ids.end()).size() != ids.size())
        throw std::invalid_argument("duplicate base id");
    for (unsigned id : ids) check_base(first, id);
    same_env(first.env_, v);
}

} // namespace

// T5.3-T5.4.
static void bind_sda_algorithms(nb::module_& m) {
    // As NDAList (T2.1): a plain list converts, list outputs take .noconvert().
    nb::bind_vector<SDAList>(m, "SDAList",
        "A list of SDA held in C++; see NDAList.");

    m.def("da_der", [](const SDA& v, unsigned id) {
        EnvGuard g(v.env_);
        check_base(v, id);
        return make_da<SDA>([&] { return da::da_der(v, id); });
    }, nb::arg("v"), nb::arg("base_id"));
    m.def("da_der", [](const SDA& v, unsigned id, SDA& result) {
        EnvGuard g(v.env_);
        same_env(v.env_, result);
        check_base(v, id);
        da::da_der(v, id, result);
    }, nb::arg("v"), nb::arg("base_id"), nb::arg("result"));
    m.def("da_int", [](const SDA& v, unsigned id) {
        EnvGuard g(v.env_);
        check_base(v, id);
        return make_da<SDA>([&] { return da::da_int(v, id); });
    }, nb::arg("v"), nb::arg("base_id"));
    m.def("da_int", [](const SDA& v, unsigned id, SDA& result) {
        EnvGuard g(v.env_);
        same_env(v.env_, result);
        check_base(v, id);
        da::da_int(v, id, result);
    }, nb::arg("v"), nb::arg("base_id"), nb::arg("result"));

    m.def("da_substitute_const", [](const SDA& iv, unsigned id, double x, SDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute_const(iv, id, x, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("x"), nb::arg("ov"));
    m.def("da_substitute", [](const SDA& iv, unsigned id, const SDA& v, SDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, v);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute(iv, id, v, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("v"), nb::arg("ov"));
    // C++ has no SDA form taking a double; it is da_substitute_const, as for NDA.
    m.def("da_substitute", [](const SDA& iv, unsigned id, double x, SDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute_const(iv, id, x, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("x"), nb::arg("ov"));
    m.def("da_substitute", [](const SDA& iv, std::vector<unsigned> ids, SDAList& v, SDA& ov) {
        EnvGuard g(iv.env_);
        check_multi(iv, ids, v);
        same_env(iv.env_, ov);
        da::da_substitute(iv, ids, v, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("v"), nb::arg("ov"));
    m.def("da_substitute", [](SDAList& ivecs, std::vector<unsigned> ids, SDAList& v,
                              SDAList& ovecs) {
        if (ivecs.size() != ovecs.size())
            throw std::invalid_argument("ivecs and ovecs must have the same length");
        if (ivecs.empty()) return;
        EnvGuard g(ivecs[0].env_);
        check_multi(ivecs[0], ids, v);
        same_env(ivecs[0].env_, ivecs);
        same_env(ivecs[0].env_, ovecs);
        check_disjoint(ovecs, ivecs);
        check_disjoint(ovecs, v);
        da::da_substitute(ivecs, ids, v, ovecs);
    }, nb::arg("ivecs"), nb::arg("base_id"), nb::arg("v"), nb::arg("ovecs").noconvert());

    m.def("da_composition", [](SDAList& ivecs, SDAList& v, SDAList& ovecs) {
        if (ivecs.size() != ovecs.size())
            throw std::invalid_argument("ivecs and ovecs must have the same length");
        if (ivecs.empty()) return;
        EnvGuard g(ivecs[0].env_);
        if (v.size() != ivecs[0].env_->layout().num_vars())
            throw std::invalid_argument("v must have nvars() entries");
        same_env(ivecs[0].env_, ivecs);
        same_env(ivecs[0].env_, v);
        same_env(ivecs[0].env_, ovecs);
        check_disjoint(ovecs, ivecs);
        check_disjoint(ovecs, v);
        da::da_composition(ivecs, v, ovecs);
    }, nb::arg("ivecs"), nb::arg("v"), nb::arg("ovecs").noconvert());
    // Returns the values (list[Expr]), as the NDA form with a float point does.
    m.def("da_composition", [](SDAList& ivecs, std::vector<double> v) {
        std::vector<Expression> out;
        if (ivecs.empty()) return out;
        EnvGuard g(ivecs[0].env_);
        if (v.size() != ivecs[0].env_->layout().num_vars())
            throw std::invalid_argument("v must have nvars() entries");
        same_env(ivecs[0].env_, ivecs);
        da::da_composition(ivecs, v, out);
        return out;
    }, nb::arg("ivecs"), nb::arg("v"));

    // Adding overloads replaces the function objects, so the T2.2 aliases
    // must be set again.
    m.attr("der") = m.attr("da_der");
    m.attr("int_") = m.attr("da_int");
    m.attr("substitute") = m.attr("da_substitute");
    m.attr("compose") = m.attr("da_composition");

    se_bridge::bind_evaluate<SDA>(m);
}

nb::class_<SDA> bind_sda(nb::module_& m, nb::class_<NDA>& nda) {
    nb::class_<SDA> c(m, "SDA",
        "A symbolic DA vector: a truncated power series whose coefficients are\n"
        "symbolic expressions (Expr). Mixed operations with NDA give an SDA.");
    c.def("__init__", [](SDA* self, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) SDA();
    }, nb::kw_only(), nb::arg("env") = nb::none(),
       "SDA(): zero. SDA(x): the constant x (float, int as an exact integer,\n"
       "Expr or symengine.Basic). In the current env, or env.");
    c.def("__init__", [](SDA* self, double x, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) SDA(x);
    }, nb::arg("x"), nb::kw_only(), nb::arg("env") = nb::none());
    c.def("__init__", [](SDA* self, nb::handle x, da::DAEnv* env) {
        Expression e = to_expr(x);
        EnvGuard g(env_or_current(env));
        new (self) SDA(e);
    }, nb::arg("x"), nb::kw_only(), nb::arg("env") = nb::none());

    c.def("copy", [](const SDA& v) { EnvGuard g(v.env_); return SDA(v); }, "A copy.");
    c.def("__copy__", [](const SDA& v) { EnvGuard g(v.env_); return SDA(v); }, "A copy.");
    c.def("__deepcopy__", [](const SDA& v, nb::handle) { EnvGuard g(v.env_); return SDA(v); },
          nb::arg("memo"), "A copy.");

    c.def_prop_rw("con",
        [](const SDA& v) { EnvGuard g(v.env_); return v.con(); },
        [](SDA& v, nb::handle x) {
            Expression e = to_expr(x);
            EnvGuard g(v.env_);
            v.reset_const_expr(e);
        },
        "The constant part, an Expr. Setting it (Expr, int, float or\n"
        "symengine.Basic) resets the vector to that constant.");
    c.def("con_symengine", [](const SDA& v) {
        Expression e = [&] { EnvGuard g(v.env_); return v.con(); }();
        return se_bridge::to_symengine(e);
    }, "The constant part as a symengine.py object.");
    c.def_prop_ro("length", [](const SDA& v) { EnvGuard g(v.env_); return v.length(); },
                  "Number of stored coefficients (up to the last non-zero one).");
    c.def_prop_ro("n_element", [](const SDA& v) { EnvGuard g(v.env_); return v.n_element(); },
                  "Number of non-zero coefficients.");
    c.def_prop_ro("nvars", [](const SDA& v) {
        EnvGuard g(v.env_);
        return v.env_->layout().num_vars();
    }, "Number of variables of the vector's env.");
    c.def_prop_ro("order", [](const SDA& v) {
        EnvGuard g(v.env_);
        return v.env_->layout().max_order();
    }, "Current truncation order of the vector's env.");
    c.def("iszero", [](const SDA& v) { EnvGuard g(v.env_); return v.iszero(); },
          "True if every coefficient is zero.");
    c.def("clean", [](SDA& v) { EnvGuard g(v.env_); v.clean(); },
          "Set the coefficients for which Expr.is_zero() holds to exactly zero.");
    c.def("reset", [](SDA& v) { EnvGuard g(v.env_); v.reset(); }, "Set every coefficient to zero.");

    c.def("element", [](SDA& v, int i) { EnvGuard g(v.env_); return v.element_expr(i); },
          nb::arg("i"),
          "element(i): coefficient i in monomial order (see exponents()), 0 past\n"
          "length. element(exps): the coefficient of the monomial with the given\n"
          "exponent of each variable. An Expr.");
    c.def("element", [](SDA& v, const std::vector<int>& exps) {
        EnvGuard g(v.env_);
        return v.element_expr(exponents(exps));
    }, nb::arg("exps"));
    c.def("set_element", [](SDA& v, const std::vector<int>& exps, nb::handle value) {
        Expression e = to_expr(value);
        EnvGuard g(v.env_);
        v.set_element(exponents(exps), e);
    }, nb::arg("exps"), nb::arg("value"),
       "Set the coefficient of the monomial with exponents exps (Expr, int,\n"
       "float or symengine.Basic).");
    c.def("index_element", [](const SDA& v, unsigned i) {
        EnvGuard g(v.env_);
        da::Layout& l = v.env_->layout();
        std::vector<unsigned> exps(l.num_vars());
        Expression value;
        da::detail::ad_elem(l, v.env_->pool<Expression>(), v.slot_, i + 1, exps.data(), value);
        return std::make_pair(std::vector<int>(exps.begin(), exps.end()), value);
    }, nb::arg("i"), "(exponents, coefficient) of coefficient i in monomial order.");
    c.def("coeffs", [](const SDA& v, bool as_symengine) {
        std::vector<Expression> cs;
        {
            EnvGuard g(v.env_);
            const Expression* p = v.env_->pool<Expression>().slot(v.slot_);
            cs.assign(p, p + v.length());
        }
        nb::typed<nb::list, nb::any> out;
        for (const Expression& e : cs)
            out.append(as_symengine ? se_bridge::to_symengine(e) : nb::cast(e));
        return out;
    }, nb::arg("as_symengine") = false,
       "The length stored coefficients: a list of Expr, or of symengine.py\n"
       "objects with as_symengine=True.");

    c.def("simplify", [](const SDA& v) {
        return map_coeffs(v, [](Expression& e) { if (!da::is_zero(e)) da::simplified_expr(e); });
    }, "A new SDA with every coefficient simplified.");
    c.def("expand", [](const SDA& v) {
        return map_coeffs(v, [](Expression& e) { e = SymEngine::expand(e); });
    }, "A new SDA with every coefficient expanded.");
    c.def("subs", [](const SDA& v, const nb::typed<nb::dict, nb::any, nb::any>& d) {
        SymEngine::map_basic_basic map = to_map(d);
        return map_coeffs(v, [&](Expression& e) { e = e.subs(map); });
    }, nb::arg("mapping"),
       "A new SDA with mapping ({symbol: value}, Expr, int, float or\n"
       "symengine.Basic) substituted in every coefficient.");

    c.def("__repr__", [](const SDA& v) {
        EnvGuard g(v.env_);
        const da::Layout& l = v.env_->layout();
        std::ostringstream os;
        os << "SDA(order=" << l.max_order() << ", nvars=" << l.num_vars()
           << ", nonzero=" << v.n_element() << ")";
        return os.str();
    }, "SDA(order=..., nvars=..., nonzero=...).");
    c.def("__str__", [](const SDA& v) {
        EnvGuard g(v.env_);
        std::ostringstream os;
        os << v;
        return os.str();
    }, "The table of terms that C++ prints.");

    // T5.2, operand order SDA, float, int, Expr, NDA, complex (T6.2),
    // symengine.Basic (A.5 rule 1). An int becomes an exact Expr integer, not
    // a float. C++ has no compound SDA-NDA operators, so those promote the NDA
    // first. SDA (op) complex gives a CSDA; it comes before symengine.Basic so
    // that a type checker, to which an untyped symengine.Basic is Any, still
    // infers CSDA for it.
    bind_arith(c);
    bind_env_property(c);
    exprops<nb::int_>(c, int_operand);
    detail_arith::binops<Expression>(c);
    detail_arith::inplace<Expression>(c, "__iadd__", [](SDA& a, const Expression& x) { a += x; });
    detail_arith::inplace<Expression>(c, "__isub__", [](SDA& a, const Expression& x) { a -= x; });
    detail_arith::inplace<Expression>(c, "__imul__", [](SDA& a, const Expression& x) { a *= x; });
    detail_arith::inplace<Expression>(c, "__itruediv__",
                                      [](SDA& a, const Expression& x) { a /= x; });
    detail_arith::binops<NDA>(c);
    detail_arith::inplace<NDA>(c, "__iadd__", [](SDA& a, const NDA& b) { a += da::promote(b); });
    detail_arith::inplace<NDA>(c, "__isub__", [](SDA& a, const NDA& b) { a -= da::promote(b); });
    detail_arith::inplace<NDA>(c, "__imul__", [](SDA& a, const NDA& b) { a *= da::promote(b); });
    detail_arith::inplace<NDA>(c, "__itruediv__",
                               [](SDA& a, const NDA& b) { a /= da::promote(b); });
    detail_arith::binops<std::complex<double>>(c);
    exprops<se_bridge::SymengineBasic>(c, symengine_operand);

    // NDA (op) Expr gives an SDA (interop.h).
    detail_arith::binops<Expression>(nda);

    m.def("promote", [](const NDA& v) {
        EnvGuard g(v.env_);
        return make_da<SDA>([&] { return da::promote(v); });
    }, nb::arg("nda"),
       "The NDA as an SDA (a CNDA as a CSDA) in the same env; integer-valued\n"
       "coefficients become exact integers.");
    m.def("svar", [](int i, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        if (i < 0 || i >= static_cast<int>(da::da_current_env().layout().num_vars()))
            throw std::out_of_range("svar: index out of range");
        return make_da<SDA>([&] { return da::promote(da::da_base(static_cast<unsigned>(i))); });
    }, nb::arg("i"), nb::kw_only(), nb::arg("env") = nb::none(),
       "promote(var(i)): the symbolic base vector of variable i.");

    m.def("sqrt", &unary<da::sqrt>);
    m.def("exp", &unary<da::exp>);
    m.def("log", &unary<da::log>);
    m.def("sin", &unary<da::sin>);
    m.def("cos", &unary<da::cos>);
    m.def("tan", &unary<da::tan>);
    m.def("asin", &unary<da::asin>);
    m.def("acos", &unary<da::acos>);
    m.def("atan", &unary<da::atan>);
    m.def("sinh", &unary<da::sinh>);
    m.def("cosh", &unary<da::cosh>);
    m.def("tanh", &unary<da::tanh>);
    m.def("erf", &unary<da::erf>);
    m.def("pow", [](const SDA& v, int n) {
        EnvGuard g(v.env_);
        return make_da<SDA>([&] { return da::pow(v, n); });
    });
    m.def("pow", [](const SDA& v, double x) {
        EnvGuard g(v.env_);
        return make_da<SDA>([&] { return da::pow(v, x); });
    });

    bind_sda_algorithms(m);
    return c;
}
