// bind_cnda.cpp — CNDA (std::complex<NDA>), CNDAList, the CNDA math functions
// and the complex algorithms (plan T3.1-T3.4).
#include "arith.h"

#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>

#include <sstream>

using Complex = std::complex<double>;
using da::get_real;
using da::get_imag;

void bind_cnda(nb::module_& m, nb::class_<NDA>& nda) {
    nb::class_<CNDA> c(m, "CNDA",
        "A complex numeric DA vector: a real and an imaginary NDA.");
    c.def("__init__", [](CNDA* self, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) CNDA();
    }, nb::kw_only(), nb::arg("env") = nb::none(),
       "CNDA(): zero. CNDA(re, im=None): re + i*im from NDA parts of one env.\n"
       "CNDA(z): the constant z. Without parts, in the current env or env.");
    c.def("__init__", [](CNDA* self, const NDA& re, const NDA* im) {
        EnvGuard g(re.env_);
        if (im) same_env(re.env_, *im);
        new (self) CNDA(re, im ? *im : NDA(0.0));
    }, nb::arg("re"), nb::arg("im") = nb::none());
    c.def("__init__", [](CNDA* self, Complex z, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) CNDA(NDA(z.real()), NDA(z.imag()));
    }, nb::arg("z"), nb::kw_only(), nb::arg("env") = nb::none());

    c.def("copy", [](const CNDA& v) {
        EnvGuard g(env_of(v));
        return make_da<CNDA>([&] { return v; });
    }, "A copy.");
    c.def("__copy__", [](const CNDA& v) {
        EnvGuard g(env_of(v));
        return make_da<CNDA>([&] { return v; });
    }, "A copy.");
    c.def("__deepcopy__", [](const CNDA& v, nb::handle) {
        EnvGuard g(env_of(v));
        return make_da<CNDA>([&] { return v; });
    }, nb::arg("memo"), "A copy.");

    c.def_prop_rw("real",
        [](const CNDA& v) {
            EnvGuard g(env_of(v));
            return make_da<NDA>([&] { return get_real(v); });
        },
        [](CNDA& v, const NDA& x) {
            EnvGuard g(env_of(v));
            same_env(env_of(v), x);
            get_real(v) = x;
        }, "The real part (a copy); assigning sets it.");
    c.def_prop_rw("imag",
        [](const CNDA& v) {
            EnvGuard g(env_of(v));
            return make_da<NDA>([&] { return get_imag(v); });
        },
        [](CNDA& v, const NDA& x) {
            EnvGuard g(env_of(v));
            same_env(env_of(v), x);
            get_imag(v) = x;
        }, "The imaginary part (a copy); assigning sets it.");

    c.def("__repr__", [](const CNDA& v) {
        EnvGuard g(env_of(v));
        const da::Layout& l = env_of(v)->layout();
        std::ostringstream os;
        os << "CNDA(order=" << l.max_order() << ", nvars=" << l.num_vars()
           << ", nonzero=(" << get_real(v).n_element() << ", " << get_imag(v).n_element()
           << "))";
        return os.str();
    }, "CNDA(order=..., nvars=..., nonzero=(re, im)).");
    c.def("__str__", [](const CNDA& v) {
        EnvGuard g(env_of(v));
        std::ostringstream os;
        os << v;
        return os.str();
    }, "The table of terms that C++ prints.");

    // Operators (arith.h), plus NDA (op) complex and complex (op) NDA from
    // da.h, which return a CNDA.
    bind_complex_arith(c);
    bind_env_property(c);
    detail_arith::binops<NDA>(c);        // CNDA (op) NDA, NDA (op) CNDA
    detail_arith::binops<Complex>(nda);

    // Like NDAList (T2.1): a plain list converts, except for output arguments.
    nb::bind_vector<CNDAList>(m, "CNDAList",
        "A list of CNDA held in C++; see NDAList.");
    bind_complex_functions<double>(m);

    m.def("read_cd_from_file", [](const std::string& filename, CNDA& cd) {
        EnvGuard g(env_of(cd));
        return da::read_cd_from_file(filename, cd);
    }, nb::arg("filename"), nb::arg("cd"),
       "Read a CNDA written by str()/print into cd; False if the file is\n"
       "incomplete.");
    m.def("compare_cd_vectors", [](CNDA& a, CNDA& b, double eps) {
        EnvGuard g(env_of(a));
        same_env(env_of(a), b);
        return da::compare_cd_vectors(a, b, eps);
    }, nb::arg("a"), nb::arg("b"), nb::arg("eps") = 1e-15,
       "compare_da_vectors of the real parts and of the imaginary parts.");
    m.def("compare_cd_with_file", [](const std::string& filename, CNDA& d, double eps) {
        EnvGuard g(env_of(d));
        return da::compare_cd_with_file(filename, d, eps);
    }, nb::arg("filename"), nb::arg("d"), nb::arg("eps") = 1e-15,
       "compare_cd_vectors of the CNDA stored in filename and d.");
}
