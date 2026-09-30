// arith.h — operator and function sets shared by the DA classes (plan T1.4, A.5).
#pragma once

#include "common.h"

#include <nanobind/stl/complex.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>

// Returns f()'s result as a new Python instance, constructed in place: the
// prvalue initializes the instance directly (guaranteed elision). Returning
// by value would move it instead, and the move ctor gives the source a fresh
// slot, which is then freed and zeroed at full length (plan T2.6 step 3).
// Typed so that the signature (and the stub) names DA, not object.
template <class DA, class F>
nb::typed<nb::object, DA> make_da(F&& f) {
    nb::object o = nb::inst_alloc(nb::type<DA>());
    new (nb::inst_ptr<DA>(o)) DA(f());
    nb::inst_mark_ready(o);
    return o;
}

// The docstring of an arithmetic operator, from its name.
inline const char* op_doc(const char* name) {
    static const char* const docs[][2] = {
        {"__add__", "Return self + value."},       {"__radd__", "Return value + self."},
        {"__sub__", "Return self - value."},       {"__rsub__", "Return value - self."},
        {"__mul__", "Return self * value."},       {"__rmul__", "Return value * self."},
        {"__truediv__", "Return self / value."},   {"__rtruediv__", "Return value / self."},
        {"__pow__", "Return self ** value."},      {"__rpow__", "Return value ** self."},
        {"__iadd__", "self += value, in place."},  {"__isub__", "self -= value, in place."},
        {"__imul__", "self *= value, in place."},  {"__itruediv__", "self /= value, in place."},
        {"__neg__", "Return -self."},              {"__pos__", "Return a copy of self."},
    };
    for (const auto& d : docs)
        if (std::strcmp(d[0], name) == 0) return d[1];
    return "";
}

namespace detail_arith {

// Returns op(a, b) as a new Python instance of op's result type.
template <class Op, class A, class B>
auto apply(const Op& op, const A& a, const B& b) {
    return make_da<std::decay_t<decltype(op(a, b))>>([&] { return op(a, b); });
}

// Binds DA (op) S and, unless S is DA, S (op) DA.
template <class S, class DA, class Op>
void binop(nb::class_<DA>& c, const char* name, const char* rname, Op op) {
    c.def(name, [op](const DA& a, const S& x) {
        EnvGuard g(env_of(a));
        return apply(op, a, x);
    }, nb::is_operator(), op_doc(name));
    if constexpr (!std::is_same_v<S, DA>)
        c.def(rname, [op](const DA& a, const S& x) {
            EnvGuard g(env_of(a));
            return apply(op, x, a);
        }, nb::is_operator(), op_doc(rname));
}

// Binds the in-place form for S; returns the same Python object.
template <class S, class DA, class Op>
void inplace(nb::class_<DA>& c, const char* name, Op op) {
    c.def(name, [op](DA& a, const S& x) -> DA& { EnvGuard g(env_of(a)); op(a, x); return a; },
          nb::is_operator(), nb::rv_policy::none, op_doc(name));
}

constexpr auto add = [](const auto& a, const auto& b) { return a + b; };
constexpr auto sub = [](const auto& a, const auto& b) { return a - b; };
constexpr auto mul = [](const auto& a, const auto& b) { return a * b; };
constexpr auto divide = [](const auto& a, const auto& b) { return a / b; };

// The four binary operators of DA with the scalar S, both sides.
template <class S, class DA>
void binops(nb::class_<DA>& c) {
    binop<S>(c, "__add__", "__radd__", add);
    binop<S>(c, "__sub__", "__rsub__", sub);
    binop<S>(c, "__mul__", "__rmul__", mul);
    binop<S>(c, "__truediv__", "__rtruediv__", divide);
}

template <class DA>
void unops(nb::class_<DA>& c) {
    c.def("__neg__", [](const DA& a) {
        EnvGuard g(env_of(a));
        return make_da<DA>([&] { return -a; });
    }, nb::is_operator(), op_doc("__neg__"));
    c.def("__pos__", [](const DA& a) {
        EnvGuard g(env_of(a));
        return make_da<DA>([&] { return +a; });
    }, nb::is_operator(), op_doc("__pos__"));
    c.def("__pow__", [](const DA& a, int n) {
        EnvGuard g(env_of(a));
        return make_da<DA>([&] { return da::pow(a, n); });
    }, nb::is_operator(), op_doc("__pow__"));
    c.def("__pow__", [](const DA& a, double x) {
        EnvGuard g(env_of(a));
        return make_da<DA>([&] { return da::pow(a, x); });
    }, nb::is_operator(), op_doc("__pow__"));
}

} // namespace detail_arith

template <class DA>
void bind_arith(nb::class_<DA>& c) {
    using namespace detail_arith;
    binops<DA>(c);
    binops<double>(c);

    inplace<DA>(c, "__iadd__", [](DA& a, const DA& b) { a += b; });
    inplace<double>(c, "__iadd__", [](DA& a, double b) { a += b; });
    inplace<DA>(c, "__isub__", [](DA& a, const DA& b) { a -= b; });
    inplace<double>(c, "__isub__", [](DA& a, double b) { a -= b; });
    inplace<DA>(c, "__imul__", [](DA& a, const DA& b) { a *= b; });
    inplace<double>(c, "__imul__", [](DA& a, double b) { a *= b; });
    inplace<DA>(c, "__itruediv__", [](DA& a, const DA& b) { a /= b; });
    inplace<double>(c, "__itruediv__", [](DA& a, double b) { a /= b; });

    unops(c);
}

// std::complex<DAVector<T>> (CNDA, CSDA) as da.h defines its operators:
// C (op) C, C (op) complex, C (op) float, both sides (no C (op) DAVector<T>).
// da.h has no compound operators for them; += and -= act on the parts in
// place, *= and /= assign the result of * and /.
template <class C>
void bind_complex_arith(nb::class_<C>& c) {
    using namespace detail_arith;
    using Complex = std::complex<double>;
    using da::get_real;
    using da::get_imag;
    binops<C>(c);
    binops<double>(c);
    binops<Complex>(c);

    inplace<C>(c, "__iadd__", [](C& a, const C& b) {
        get_real(a) += get_real(b);
        get_imag(a) += get_imag(b);
    });
    inplace<double>(c, "__iadd__", [](C& a, double x) { get_real(a) += x; });
    inplace<Complex>(c, "__iadd__", [](C& a, Complex z) {
        get_real(a) += z.real();
        get_imag(a) += z.imag();
    });
    inplace<C>(c, "__isub__", [](C& a, const C& b) {
        get_real(a) -= get_real(b);
        get_imag(a) -= get_imag(b);
    });
    inplace<double>(c, "__isub__", [](C& a, double x) { get_real(a) -= x; });
    inplace<Complex>(c, "__isub__", [](C& a, Complex z) {
        get_real(a) -= z.real();
        get_imag(a) -= z.imag();
    });
    inplace<C>(c, "__imul__", [](C& a, const C& b) { a = a * b; });
    inplace<double>(c, "__imul__", [](C& a, double x) { a = a * x; });
    inplace<Complex>(c, "__imul__", [](C& a, Complex z) { a = a * z; });
    inplace<C>(c, "__itruediv__", [](C& a, const C& b) { a = divide(a, b); });
    inplace<double>(c, "__itruediv__", [](C& a, double x) { a = divide(a, x); });
    inplace<Complex>(c, "__itruediv__", [](C& a, Complex z) { a = divide(a, z); });

    unops(c);
}

namespace detail_arith {

template <class C, C (*F)(const C&)>
auto complex_unary(const C& c) {
    EnvGuard g(env_of(c));
    return make_da<C>([&] { return F(c); });
}

// The env of a cd_composition call: v's, which has one entry per variable.
template <class V>
da::DAEnv* composition_env(const V& v) {
    if (v.empty()) throw std::invalid_argument("v must have nvars() entries");
    return env_of(v[0]);
}

// Checks shared by the cd_composition forms; C++ only asserts them.
template <class In, class V, class Out>
void check_composition(da::DAEnv* e, const In& ivecs, const V& v, const Out& ovecs) {
    if (v.size() != e->layout().num_vars())
        throw std::invalid_argument("v must have nvars() entries");
    if (ivecs.size() != ovecs.size())
        throw std::invalid_argument("ivecs and ovecs must have the same length");
    same_env(e, ivecs);
    same_env(e, v);
    same_env(e, ovecs);
}

} // namespace detail_arith

// The functions and cd_composition forms of std::complex<DAVector<T>>, added
// to the names after the real overloads. functions.h has no complex sin, cos,
// tan, sinh, cosh, tanh or erf. abs is a double for CNDA (the larger part
// norm) and an SDA for CSDA (sqrt(re^2 + im^2)).
template <class T>
void bind_complex_functions(nb::module_& m) {
    using namespace detail_arith;
    using D = da::DAVector<T>;
    using C = std::complex<D>;
    using DList = std::vector<D>;
    using CList = std::vector<C>;
    m.def("sqrt", &complex_unary<C, da::sqrt<T>>);
    m.def("exp", &complex_unary<C, da::exp<T>>);
    m.def("log", &complex_unary<C, da::log<T>>);
    m.def("asin", &complex_unary<C, da::asin<T>>);
    m.def("acos", &complex_unary<C, da::acos<T>>);
    m.def("atan", &complex_unary<C, da::atan<T>>);
    m.def("asinh", &complex_unary<C, da::asinh<T>>);
    m.def("acosh", &complex_unary<C, da::acosh<T>>);
    m.def("atanh", &complex_unary<C, da::atanh<T>>);
    if constexpr (std::is_same_v<T, double>)
        m.def("abs", [](const C& v) { EnvGuard g(env_of(v)); return da::abs(v); },
              "abs(CNDA): the larger of norm() of the real and the imaginary part.");
    else
        m.def("abs", [](const C& v) {
            EnvGuard g(env_of(v));
            return make_da<D>([&] { return da::abs(v); });
        }, "abs(CSDA): the modulus sqrt(re**2 + im**2), as an SDA.");
    m.def("pow", [](const C& v, int n) {
        EnvGuard g(env_of(v));
        return make_da<C>([&] { return da::pow(v, n); });
    });
    m.def("pow", [](const C& v, double x) {
        EnvGuard g(env_of(v));
        return make_da<C>([&] { return da::pow(v, x); });
    });

    // Only the real-ivecs form reads v after it has reset ovecs, so only
    // there may an output not be an input.
    m.def("cd_composition", [](DList& ivecs, CList& v, CList& ovecs) {
        da::DAEnv* e = composition_env(v);
        EnvGuard g(e);
        check_composition(e, ivecs, v, ovecs);
        for (const C& o : ovecs)
            for (const C& i : v)
                if (da::get_real(o).slot_ == da::get_real(i).slot_)
                    throw std::invalid_argument("an output vector is also an input vector");
        da::cd_composition(ivecs, v, ovecs);
    }, nb::arg("ivecs"), nb::arg("v"), nb::arg("ovecs").noconvert(),
       "Compose a map with complex DA vectors: ovecs[i] = ivecs[i](v[0], ..., v[nvars-1]).\n\n"
       "ivecs is a real or complex map, v a real or complex list with nvars()\n"
       "entries; the results go to ovecs (a CNDAList or CSDAList of len(ivecs)).");
    m.def("cd_composition", [](CList& ivecs, CList& v, CList& ovecs) {
        da::DAEnv* e = composition_env(v);
        EnvGuard g(e);
        check_composition(e, ivecs, v, ovecs);
        da::cd_composition(ivecs, v, ovecs);
    }, nb::arg("ivecs"), nb::arg("v"), nb::arg("ovecs").noconvert());
    m.def("cd_composition", [](CList& ivecs, DList& v, CList& ovecs) {
        da::DAEnv* e = composition_env(v);
        EnvGuard g(e);
        check_composition(e, ivecs, v, ovecs);
        da::cd_composition(ivecs, v, ovecs);
    }, nb::arg("ivecs"), nb::arg("v"), nb::arg("ovecs").noconvert());
}
