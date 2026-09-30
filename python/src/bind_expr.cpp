// bind_expr.cpp — Expr (SymEngine::Expression), symbols() and the symengine.py
// interop (plan T4.1, T4.3, T4.4).
#include "se_bridge.h"

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "da/symbolic_ops.h"

#include <symengine/visitor.h>

#include <sstream>

using SymEngine::Expression;
using se_bridge::from_int;

namespace {

Expression parse(const std::string& s) {
    try {
        return Expression(SymEngine::parse(s));
    } catch (const SymEngine::ParseError& e) {
        throw std::invalid_argument(std::string("cannot parse Expr: ") + e.what());
    }
}

std::string to_string(const Expression& e) {
    std::ostringstream os;
    os << e;
    return os.str();
}

// Expr, int or float (subs values, dict keys).
Expression to_expr(nb::handle h) {
    if (nb::isinstance<Expression>(h)) return nb::cast<const Expression&>(h);
    if (PyLong_Check(h.ptr())) return from_int(nb::borrow<nb::int_>(h));
    if (PyFloat_Check(h.ptr())) return Expression(PyFloat_AsDouble(h.ptr()));
    throw nb::type_error("expected Expr, int or float");
}

// Binds Expr (op) S and S (op) Expr for S = Expr, int, float, in that order.
template <class Op>
void binop(nb::class_<Expression>& c, const char* name, const char* rname, Op op) {
    c.def(name, [op](const Expression& a, const Expression& b) { return op(a, b); },
          nb::is_operator(), op_doc(name));
    c.def(name, [op](const Expression& a, const nb::int_& b) { return op(a, from_int(b)); },
          nb::is_operator(), op_doc(name));
    c.def(rname, [op](const Expression& a, const nb::int_& b) { return op(from_int(b), a); },
          nb::is_operator(), op_doc(rname));
    c.def(name, [op](const Expression& a, double b) { return op(a, Expression(b)); },
          nb::is_operator(), op_doc(name));
    c.def(rname, [op](const Expression& a, double b) { return op(Expression(b), a); },
          nb::is_operator(), op_doc(rname));
}

} // namespace

void bind_expr(nb::module_& m) {
    nb::class_<Expression> c(m, "Expr",
        "A symbolic scalar (SymEngine::Expression), the coefficient type of SDA.");
    c.def("__init__", [](Expression* self, const nb::int_& i) { new (self) Expression(from_int(i)); },
          "Expr(int): an exact integer. Expr(float). Expr(str): the parsed\n"
          "expression, e.g. Expr('a + 2*b').");
    c.def("__init__", [](Expression* self, double x) { new (self) Expression(x); });
    c.def("__init__", [](Expression* self, const std::string& s) { new (self) Expression(parse(s)); });

    binop(c, "__add__", "__radd__", [](const Expression& a, const Expression& b) { return a + b; });
    binop(c, "__sub__", "__rsub__", [](const Expression& a, const Expression& b) { return a - b; });
    binop(c, "__mul__", "__rmul__", [](const Expression& a, const Expression& b) { return a * b; });
    binop(c, "__truediv__", "__rtruediv__",
          [](const Expression& a, const Expression& b) { return a / b; });
    binop(c, "__pow__", "__rpow__",
          [](const Expression& a, const Expression& b) { return SymEngine::pow(a, b); });
    c.def("__neg__", [](const Expression& a) { return -a; }, nb::is_operator(), op_doc("__neg__"));
    c.def("__pos__", [](const Expression& a) { return a; }, nb::is_operator(), op_doc("__pos__"));

    // Structural equality with an Expr, int or float; any other object is
    // unequal (the last overloads, which also keep the stub's __eq__
    // compatible with object.__eq__).
    const char* eq_doc = "Structural equality with an Expr, int or float.";
    const char* ne_doc = "Negation of __eq__.";
    c.def("__eq__", [](const Expression& a, const Expression& b) { return a == b; },
          nb::is_operator(), eq_doc);
    c.def("__eq__", [](const Expression& a, const nb::int_& b) { return a == from_int(b); },
          nb::is_operator(), eq_doc);
    c.def("__eq__", [](const Expression& a, double b) { return a == Expression(b); },
          nb::is_operator(), eq_doc);
    c.def("__eq__", [](const Expression&, nb::handle) { return false; }, nb::is_operator(), eq_doc);
    c.def("__ne__", [](const Expression& a, const Expression& b) { return !(a == b); },
          nb::is_operator(), ne_doc);
    c.def("__ne__", [](const Expression& a, const nb::int_& b) { return !(a == from_int(b)); },
          nb::is_operator(), ne_doc);
    c.def("__ne__", [](const Expression& a, double b) { return !(a == Expression(b)); },
          nb::is_operator(), ne_doc);
    c.def("__ne__", [](const Expression&, nb::handle) { return true; }, nb::is_operator(), ne_doc);
    c.def("__hash__", [](const Expression& a) { return a.get_basic()->hash(); },
          "SymEngine's hash, consistent with __eq__.");

    c.def("__str__", &to_string, "The expression as text, e.g. 'a + 2*b'.");
    c.def("__repr__", [](const Expression& a) {
        return "Expr(" + nb::cast<std::string>(nb::repr(nb::str(to_string(a).c_str()))) + ")";
    }, "Expr('...').");
    c.def("__float__", [](const Expression& a) {
        if (!SymEngine::free_symbols(*a.get_basic()).empty())
            throw nb::type_error("cannot convert an Expr with free symbols to float");
        return SymEngine::eval_double(*a.get_basic());
    }, "The numeric value; TypeError if free symbols remain.");

    c.def("subs", [](const Expression& a, const nb::typed<nb::dict, nb::any, nb::any>& d) {
        SymEngine::map_basic_basic map;
        for (auto [k, v] : d) map[to_expr(k).get_basic()] = to_expr(v).get_basic();
        return a.subs(map);
    }, nb::arg("mapping"),
       "A new Expr with mapping ({symbol: value}, Expr, int or float)\n"
       "substituted.");
    c.def("expand", [](const Expression& a) { return Expression(SymEngine::expand(a)); },
          "The expanded expression.");
    c.def("diff", [](const Expression& a, const Expression& sym) {
        if (!SymEngine::is_a<SymEngine::Symbol>(*sym.get_basic()))
            throw std::invalid_argument("diff: the variable must be a symbol");
        return a.diff(SymEngine::rcp_static_cast<const SymEngine::Symbol>(sym.get_basic()));
    }, nb::arg("sym"), "The derivative with respect to the symbol sym.");
    c.def("free_symbols", [](const Expression& a) {
        nb::typed<nb::set, Expression> out;
        for (const auto& s : SymEngine::free_symbols(*a.get_basic())) out.add(nb::cast(Expression(s)));
        return out;
    }, "The set of symbols in the expression.");
    c.def("is_zero", [](const Expression& a) { return da::is_zero(a); },
          "True if the expression is zero (the test SDA uses).");
    c.def("simplify", [](const Expression& a) {
        Expression r = a;
        da::simplified_expr(r);
        return r;
    }, "A simplified copy (the simplification SDA.simplify() uses).");

    c.def("to_symengine", [](const Expression& a) { return se_bridge::to_symengine(a); },
          "The expression as a symengine.py object: the same C++ object when\n"
          "symengine.py shares MiraDAC's libsymengine, else via a string (see\n"
          "symengine_interop_status()).");
    c.def_static("from_symengine", [](nb::handle obj) { return se_bridge::from_symengine(obj); },
                 nb::arg("obj"), "The Expr of a symengine.Basic (see to_symengine()).");

    m.def("symbols", [](const std::string& names) {
        std::vector<Expression> out;
        std::string name;
        for (char ch : names + " ") {
            if (ch == ' ' || ch == ',' || ch == '\t' || ch == '\n') {
                if (!name.empty()) out.emplace_back(SymEngine::symbol(name));
                name.clear();
            } else {
                name += ch;
            }
        }
        nb::typed<nb::tuple, Expression, nb::ellipsis> t =
            nb::steal<nb::tuple>(PyTuple_New(out.size()));
        for (std::size_t i = 0; i < out.size(); ++i)
            PyTuple_SET_ITEM(t.ptr(), i, nb::cast(out[i]).release().ptr());
        return t;
    }, nb::arg("names"),
       "Symbols named in a string separated by spaces or commas:\n"
       "a, b = symbols('a b').");

    m.def("symengine_interop_status", [] {
        const se_bridge::Status& s = se_bridge::status();
        nb::typed<nb::dict, nb::str, nb::any> d;
        d["mode"] = s.shared ? "shared" : "string";
        d["symengine_version"] = s.imported ? nb::cast(s.version) : nb::none();
        d["wrapper_path"] = s.imported ? nb::cast(s.wrapper_path) : nb::none();
        d["wrapper_needed"] = s.needed;
        d["loaded_libsymengine"] = s.loaded;
        d["expected_libsymengine"] = DA_SYMENGINE_LIB;
        d["layout_selftest"] = s.selftest;
        d["reason"] = s.reason.empty() ? nb::none() : nb::cast(s.reason);
        return d;
    }, "Whether symengine.py objects cross by pointer (mode 'shared') or by\n"
       "string (mode 'string'), with the result of each check and, in string\n"
       "mode, the first failed check (reason). The checks run once.");

    // For the import check in __init__.py and the tests.
    m.def("_libsymengine_path", &se_bridge::our_libsymengine,
          "Realpath of the libsymengine this module is linked to.");
    m.attr("_LIBSYMENGINE_EXPECTED") = DA_SYMENGINE_LIB;
    m.attr("_LIBSYMENGINE_SHA256") = DA_SYMENGINE_LIB_SHA256;
    m.def("_rcp_address", [](nb::handle obj) {
        return reinterpret_cast<std::uintptr_t>(se_bridge::rcp_for_test(obj).get());
    }, nb::arg("obj"), "Address of the SymEngine object behind obj (for tests).");
    m.def("_rcp_use_count", [](nb::handle obj) { return se_bridge::rcp_for_test(obj)->use_count(); },
          nb::arg("obj"), "SymEngine reference count of the object behind obj (for tests).");
}
