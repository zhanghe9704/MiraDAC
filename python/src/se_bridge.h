// se_bridge.h — symengine.py interop: the checks that decide between zero-copy
// ("shared") and string conversion, and the conversions (plan A.8, A.9, T4.3, T4.4).
#pragma once

#include "arith.h"

#include <symengine/expression.h>
#include <symengine/eval_double.h>
#include <symengine/parser.h>

#include <dlfcn.h>
#include <elf.h>
#include <limits.h>
#include <stdlib.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

// Set by python/CMakeLists.txt from the pin and the SymEngine found at configure time.
#ifndef DA_SYMENGINE_LIB
#error "DA_SYMENGINE_LIB must be defined"
#endif

namespace se_bridge {

using SymEngine::Expression;
using RCPBasic = SymEngine::RCP<const SymEngine::Basic>;
using C2py = PyObject* (*)(RCPBasic);

// symengine.py's `cdef class Basic` is PyObject_HEAD followed by one RCP (A.9).
static_assert(sizeof(RCPBasic) == sizeof(void*), "SymEngine RCP must be a single pointer");

// An exact integer (beyond long long through the parser).
inline Expression from_int(const nb::int_& i) {
    int overflow = 0;
    long long v = PyLong_AsLongLongAndOverflow(i.ptr(), &overflow);
    if (overflow) return Expression(SymEngine::parse(nb::cast<std::string>(nb::str(i))));
    if (v == -1 && PyErr_Occurred()) throw nb::python_error();
    return Expression(SymEngine::integer(static_cast<long>(v)));
}

inline const RCPBasic& rcp_of(PyObject* o) {
    return *reinterpret_cast<const RCPBasic*>(reinterpret_cast<const char*>(o) + sizeof(PyObject));
}

inline std::string real_path(const std::string& p) {
    char buf[PATH_MAX];
    return ::realpath(p.c_str(), buf) ? std::string(buf) : p;
}

// Realpath of the libsymengine this module is linked to.
inline std::string our_libsymengine() {
    Dl_info info{};
    auto f = static_cast<double (*)(const SymEngine::Basic&)>(&SymEngine::eval_double);
    if (!dladdr(reinterpret_cast<void*>(f), &info) || !info.dli_fname) return "";
    return real_path(info.dli_fname);
}

// DT_NEEDED entries of an ELF64 file (empty if it cannot be read).
inline std::vector<std::string> elf_needed(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::string d((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    std::vector<std::string> out;
    if (d.size() < sizeof(Elf64_Ehdr) || d.compare(0, SELFMAG, ELFMAG) != 0 ||
        d[EI_CLASS] != ELFCLASS64)
        return out;
    auto eh = reinterpret_cast<const Elf64_Ehdr*>(d.data());
    if (eh->e_shoff + std::size_t(eh->e_shnum) * sizeof(Elf64_Shdr) > d.size()) return out;
    auto sh = reinterpret_cast<const Elf64_Shdr*>(d.data() + eh->e_shoff);
    for (unsigned i = 0; i < eh->e_shnum; ++i) {
        if (sh[i].sh_type != SHT_DYNAMIC || sh[i].sh_link >= eh->e_shnum) continue;
        const Elf64_Shdr& str = sh[sh[i].sh_link];
        if (sh[i].sh_offset + sh[i].sh_size > d.size() || str.sh_offset + str.sh_size > d.size())
            return out;
        auto dyn = reinterpret_cast<const Elf64_Dyn*>(d.data() + sh[i].sh_offset);
        for (std::size_t k = 0; k < sh[i].sh_size / sizeof(Elf64_Dyn); ++k)
            if (dyn[k].d_tag == DT_NEEDED && dyn[k].d_un.d_val < str.sh_size)
                out.emplace_back(d.data() + str.sh_offset + dyn[k].d_un.d_val);
    }
    return out;
}

// Realpaths of every mapped file whose name starts with "libsymengine".
inline std::vector<std::string> mapped_libsymengine() {
    std::ifstream maps("/proc/self/maps");
    std::vector<std::string> out;
    std::string line;
    while (std::getline(maps, line)) {
        auto p = line.find('/');
        if (p == std::string::npos) continue;
        std::string path = line.substr(p);
        std::string base = path.substr(path.rfind('/') + 1);
        if (base.rfind("libsymengine", 0) != 0) continue;
        path = real_path(path);
        if (std::find(out.begin(), out.end(), path) == out.end()) out.push_back(path);
    }
    return out;
}

struct Status {
    bool shared = false;
    bool imported = false;
    std::string version, wrapper_path, reason;
    std::vector<std::string> needed, loaded;
    bool selftest = false;
    C2py c2py = nullptr;
    PyObject* basic_type = nullptr;  // symengine_wrapper.Basic, kept alive for the process
};

// Runs the four "first interop call" checks of A.8, in order.
inline Status run_checks() {
    Status s;
    nb::object w;
    try {
        nb::object se = nb::module_::import_("symengine");
        w = nb::module_::import_("symengine.lib.symengine_wrapper");
        s.imported = true;
        s.version = nb::cast<std::string>(nb::str(se.attr("__version__")));
        s.wrapper_path = real_path(nb::cast<std::string>(w.attr("__file__")));
    } catch (nb::python_error& e) {
        s.reason = std::string("symengine.py is not importable: ") + nb::str(e.value()).c_str();
        return s;
    }
    s.needed = elf_needed(s.wrapper_path);
    s.loaded = mapped_libsymengine();
    const std::string expected = DA_SYMENGINE_LIB;

    if (s.version != DA_SYMENGINE_PY_VERSION) {
        s.reason = "symengine.py is version " + s.version + ", the pin is " DA_SYMENGINE_PY_VERSION;
    } else if (std::find(s.needed.begin(), s.needed.end(), DA_SYMENGINE_SONAME) == s.needed.end()) {
        s.reason = s.wrapper_path + " does not link " DA_SYMENGINE_SONAME
                   " (it has its own SymEngine copy, as the PyPI wheel does)";
    } else if (s.loaded.size() != 1 || s.loaded[0] != expected) {
        s.reason = "the process maps " + std::to_string(s.loaded.size()) +
                   " libsymengine file(s), expected only " + expected;
    } else {
        nb::object cap = w.attr("__pyx_capi__")["c2py"];
        s.c2py = reinterpret_cast<C2py>(
            PyCapsule_GetPointer(cap.ptr(), PyCapsule_GetName(cap.ptr())));
        if (!s.c2py) throw nb::python_error();
        s.basic_type = w.attr("Basic").inc_ref().ptr();
        RCPBasic sym = SymEngine::symbol("miradac_layout_selftest");
        PyObject* o = s.c2py(sym);
        if (!o) throw nb::python_error();
        nb::object obj = nb::steal(o);
        s.selftest = PyObject_IsInstance(o, s.basic_type) == 1 &&
                     rcp_of(o).get() == sym.get() && sym->use_count() == 2;
        if (!s.selftest) s.reason = "layout self-test failed: symengine.Basic no longer holds one RCP";
    }
    s.shared = s.reason.empty();
    return s;
}

// Checks run once, on the first interop call; the result is kept for the process.
inline const Status& status() {
    static const Status* s = new Status(run_checks());
    return *s;
}

// String mode: raise with MIRADAC_REQUIRE_SHARED_SYMENGINE=1, else warn once.
inline void enter_string_mode(const Status& s) {
    const char* req = std::getenv("MIRADAC_REQUIRE_SHARED_SYMENGINE");
    const std::string msg = "miradac: symengine.py interop is in string mode (" + s.reason +
        "); build symengine.py against the pinned SymEngine with scripts/setup_symengine.sh "
        "(plan TP.2) for zero-copy conversion";
    if (req && std::strcmp(req, "1") == 0) throw std::runtime_error(msg);
    static bool warned = false;
    if (!warned) {
        warned = true;
        if (PyErr_WarnEx(PyExc_RuntimeWarning, msg.c_str(), 1) < 0) throw nb::python_error();
    }
}

inline nb::object to_symengine(const Expression& e) {
    const Status& s = status();
    if (s.shared) {
        PyObject* o = s.c2py(e.get_basic());
        if (!o) throw nb::python_error();
        return nb::steal(o);
    }
    enter_string_mode(s);
    std::ostringstream os;
    os << e;
    return nb::module_::import_("symengine").attr("sympify")(os.str());
}

inline Expression from_symengine(nb::handle obj) {
    const Status& s = status();
    if (s.shared) {
        if (PyObject_IsInstance(obj.ptr(), s.basic_type) != 1)
            throw nb::type_error("expected a symengine.Basic");
        return Expression(rcp_of(obj.ptr()));
    }
    enter_string_mode(s);
    nb::object basic = nb::module_::import_("symengine").attr("Basic");
    if (!nb::isinstance(obj, basic)) throw nb::type_error("expected a symengine.Basic");
    return Expression(SymEngine::parse(nb::cast<std::string>(nb::str(obj))));
}

// True if h is a symengine.Basic. Does not import symengine.py or run the
// checks: if symengine.py was never imported, h cannot be one of its objects.
inline bool is_symengine(nb::handle h) {
    nb::object w = nb::steal(PyImport_GetModule(nb::str("symengine.lib.symengine_wrapper").ptr()));
    if (!w.is_valid()) {
        PyErr_Clear();
        return false;
    }
    return nb::isinstance(h, w.attr("Basic"));
}

// A symengine.Basic argument, named so in signatures and stubs. Any other
// object fails the check, so overload resolution moves on.
class SymengineBasic : public nb::object {
    static bool check(PyObject* o) { return is_symengine(o); }

public:
    NB_OBJECT_DEFAULT(SymengineBasic, nb::object, "symengine.Basic", check)
};

// The SymEngine object behind an Expr or (shared mode only) a symengine.Basic.
inline const RCPBasic& rcp_for_test(nb::handle obj) {
    if (nb::isinstance<Expression>(obj)) return nb::cast<const Expression&>(obj).get_basic();
    const Status& s = status();
    if (s.shared && PyObject_IsInstance(obj.ptr(), s.basic_type) == 1) return rcp_of(obj.ptr());
    throw nb::type_error("expected an Expr, or a symengine.Basic in shared mode");
}

// Expr, int, float or symengine.Basic (T5.5).
inline Expression to_expr(nb::handle h) {
    if (nb::isinstance<Expression>(h)) return nb::cast<const Expression&>(h);
    if (PyLong_Check(h.ptr())) return from_int(nb::borrow<nb::int_>(h));
    if (PyFloat_Check(h.ptr())) return Expression(PyFloat_AsDouble(h.ptr()));
    if (is_symengine(h)) return from_symengine(h);
    throw nb::type_error("expected Expr, int, float or symengine.Basic");
}

// da::evaluate of an SDA (giving an NDA) or a CSDA (giving a CNDA).
template <class S>
auto evaluate_at(const S& s, const std::vector<RCPBasic>& syms,
                       const std::vector<double>& vals) {
    if (syms.size() != vals.size())
        throw std::invalid_argument("evaluate: syms and vals must have the same length");
    EnvGuard g(env_of(s));
    try {
        return make_da<std::decay_t<decltype(da::evaluate(s, syms, vals))>>(
            [&] { return da::evaluate(s, syms, vals); });
    } catch (const SymEngine::SymEngineException& e) {
        throw std::invalid_argument(std::string("evaluate: ") + e.what());
    }
}

// Both forms use the vector overload of da::evaluate (interop.h), which is
// faster than the map one; a dict is split into keys and values.
template <class S>
void bind_evaluate(nb::module_& m) {
    m.def("evaluate", [](const S& s, const nb::typed<nb::dict, nb::any, double>& values) {
        std::vector<RCPBasic> syms;
        std::vector<double> vals;
        for (auto [k, v] : values) {
            syms.push_back(to_expr(k).get_basic());
            vals.push_back(nb::cast<double>(v));
        }
        return evaluate_at(s, syms, vals);
    }, nb::arg("sda"), nb::arg("values"),
       "Substitute numbers for the symbols: an SDA gives an NDA, a CSDA a CNDA.\n\n"
       "values is a dict {symbol: float}, or give syms and vals as two\n"
       "sequences. Symbols are Expr or symengine.Basic; every symbol in a\n"
       "coefficient needs a value (ValueError otherwise).");
    m.def("evaluate", [](const S& s, const nb::typed<nb::sequence, nb::any>& syms,
                         std::vector<double> vals) {
        std::vector<RCPBasic> rcps;
        for (nb::handle h : syms) rcps.push_back(to_expr(h).get_basic());
        return evaluate_at(s, rcps, vals);
    }, nb::arg("sda"), nb::arg("syms"), nb::arg("vals"));
}

} // namespace se_bridge
