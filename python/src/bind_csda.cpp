// bind_csda.cpp — CSDA (std::complex<SDA>), CSDAList, promote(CNDA), the CSDA
// math functions, cd_composition and evaluate() (plan T6.1-T6.3).
//
// Operators, both sides: CSDA (op) CSDA, complex, float (da.h's is_da_coeff
// templates), SDA (std::complex's own templates), Expr and symengine.Basic
// (da.h). SDA (op) complex, which gives a CSDA, is bound in bind_sda.cpp.
#include "arith.h"
#include "se_bridge.h"

#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <sstream>

using Complex = std::complex<double>;
using da::get_real;
using da::get_imag;

void bind_csda(nb::module_& m) {
    nb::class_<CSDA> c(m, "CSDA",
        "A complex symbolic DA vector: a real and an imaginary SDA.");
    c.def("__init__", [](CSDA* self, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) CSDA();
    }, nb::kw_only(), nb::arg("env") = nb::none(),
       "CSDA(): zero, in the current env or env. CSDA(re, im=None): re + i*im\n"
       "from SDA parts of one env.");
    c.def("__init__", [](CSDA* self, const SDA& re, const SDA* im) {
        EnvGuard g(re.env_);
        if (im) same_env(re.env_, *im);
        new (self) CSDA(re, im ? *im : SDA());
    }, nb::arg("re"), nb::arg("im") = nb::none());

    c.def("copy", [](const CSDA& v) {
        EnvGuard g(env_of(v));
        return make_da<CSDA>([&] { return v; });
    }, "A copy.");
    c.def("__copy__", [](const CSDA& v) {
        EnvGuard g(env_of(v));
        return make_da<CSDA>([&] { return v; });
    }, "A copy.");
    c.def("__deepcopy__", [](const CSDA& v, nb::handle) {
        EnvGuard g(env_of(v));
        return make_da<CSDA>([&] { return v; });
    }, nb::arg("memo"), "A copy.");

    c.def_prop_rw("real",
        [](const CSDA& v) {
            EnvGuard g(env_of(v));
            return make_da<SDA>([&] { return get_real(v); });
        },
        [](CSDA& v, const SDA& x) {
            EnvGuard g(env_of(v));
            same_env(env_of(v), x);
            get_real(v) = x;
        }, "The real part (a copy); assigning sets it.");
    c.def_prop_rw("imag",
        [](const CSDA& v) {
            EnvGuard g(env_of(v));
            return make_da<SDA>([&] { return get_imag(v); });
        },
        [](CSDA& v, const SDA& x) {
            EnvGuard g(env_of(v));
            same_env(env_of(v), x);
            get_imag(v) = x;
        }, "The imaginary part (a copy); assigning sets it.");

    // No __str__: C++ has no operator<< for std::complex<SDA> (print .real
    // and .imag).
    c.def("__repr__", [](const CSDA& v) {
        EnvGuard g(env_of(v));
        const da::Layout& l = env_of(v)->layout();
        std::ostringstream os;
        os << "CSDA(order=" << l.max_order() << ", nvars=" << l.num_vars()
           << ", nonzero=(" << get_real(v).n_element() << ", " << get_imag(v).n_element()
           << "))";
        return os.str();
    }, "CSDA(order=..., nvars=..., nonzero=(re, im)).");

    bind_complex_arith(c);
    bind_env_property(c);
    detail_arith::binops<SDA>(c);
    detail_arith::binops<SymEngine::Expression>(c);
    {
        using namespace detail_arith;
        auto se_op = [&](const char* name, const char* rname, auto op) {
            c.def(name, [op](const CSDA& a, const se_bridge::SymengineBasic& x) {
                SymEngine::Expression e = se_bridge::from_symengine(x);
                EnvGuard g(env_of(a));
                return make_da<CSDA>([&] { return op(a, e); });
            }, nb::is_operator(), op_doc(name));
            c.def(rname, [op](const CSDA& a, const se_bridge::SymengineBasic& x) {
                SymEngine::Expression e = se_bridge::from_symengine(x);
                EnvGuard g(env_of(a));
                return make_da<CSDA>([&] { return op(e, a); });
            }, nb::is_operator(), op_doc(rname));
        };
        se_op("__add__", "__radd__", add);
        se_op("__sub__", "__rsub__", sub);
        se_op("__mul__", "__rmul__", mul);
        se_op("__truediv__", "__rtruediv__", divide);
    }
    c.def("__str__", [](const CSDA& v) {
        EnvGuard g(env_of(v));
        std::ostringstream os;
        os << v;
        return os.str();
    }, "The real part's table of terms, then the imaginary part's.");

    m.def("promote", [](const CNDA& v) {
        EnvGuard g(env_of(v));
        return make_da<CSDA>([&] { return da::promote(v); });
    }, nb::arg("cnda"));

    // As SDAList (T5.3): a plain list converts, list outputs take .noconvert().
    nb::bind_vector<CSDAList>(m, "CSDAList",
        "A list of CSDA held in C++; see NDAList.");
    bind_complex_functions<SymEngine::Expression>(m);
    se_bridge::bind_evaluate<CSDA>(m);
}
