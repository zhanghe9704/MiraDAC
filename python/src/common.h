// common.h — shared helpers for the miradac._core bindings (plan A.4, A.10).
#pragma once

#include <nanobind/nanobind.h>

#include "da/da.h"

#include <complex>
#include <stdexcept>
#include <vector>

namespace nb = nanobind;

using da::NDA;
using da::CNDA;
using NDAList = std::vector<NDA>;
using CNDAList = std::vector<CNDA>;

// Bound as classes. For CNDA (and CSDA) this also overrides the std::complex
// caster of nanobind/stl/complex.h, which would otherwise claim std::complex<NDA>.
NB_MAKE_OPAQUE(NDAList)
NB_MAKE_OPAQUE(CNDA)
NB_MAKE_OPAQUE(CNDAList)
#ifdef DA_WITH_SYMBOLIC
using da::SDA;
using da::CSDA;
using SDAList = std::vector<SDA>;
using CSDAList = std::vector<CSDA>;
NB_MAKE_OPAQUE(SDAList)
NB_MAKE_OPAQUE(CSDA)
NB_MAKE_OPAQUE(CSDAList)
#endif

struct EnvError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Selects a DA object's env for the duration of a bound call.
struct EnvGuard {
    da::DAEnv* prev;
    explicit EnvGuard(da::DAEnv* e) {
        if (e->retired()) throw EnvError("DA environment has been cleared");
        prev = da::da_exchange_env(e);
    }
    ~EnvGuard() { da::da_exchange_env(prev); }
    EnvGuard(const EnvGuard&) = delete;
    EnvGuard& operator=(const EnvGuard&) = delete;
};

inline da::DAEnv* current_env_or_throw() {
    da::DAEnv* e = da::da_exchange_env(nullptr);
    da::da_exchange_env(e);
    if (!e) throw EnvError("no DA environment: call miradac.init() first");
    if (e->retired()) throw EnvError("DA environment has been cleared");
    return e;
}

// The env of an env= keyword (plan T7.3), else the current env.
inline da::DAEnv* env_or_current(da::DAEnv* e) { return e ? e : current_env_or_throw(); }

template <class T>
da::DAEnv* env_of(const da::DAVector<T>& v) { return v.env_; }
template <class T>
da::DAEnv* env_of(const std::complex<da::DAVector<T>>& c) { return da::get_real(c).env_; }

// The .env property: the owning Env (bind_env.cpp never lets an env be deleted).
template <class DA>
void bind_env_property(nb::class_<DA>& c) {
    c.def_prop_ro("env", [](const DA& v) { return env_of(v); }, nb::rv_policy::reference);
}

// The engine indexes one pool with the slots of every vector it is given,
// so all of them must live in env e (C++ checks this only for operators).
template <class T>
void same_env(da::DAEnv* e, const da::DAVector<T>& v) {
    if (v.env_ != e) throw EnvError("DA vectors belong to different environments");
}
template <class T>
void same_env(da::DAEnv* e, const std::vector<da::DAVector<T>>& l) {
    for (const da::DAVector<T>& v : l) same_env(e, v);
}
template <class T>
void same_env(da::DAEnv* e, const std::complex<da::DAVector<T>>& c) {
    same_env(e, da::get_real(c));
    same_env(e, da::get_imag(c));
}
template <class T>
void same_env(da::DAEnv* e, const std::vector<std::complex<da::DAVector<T>>>& l) {
    for (const std::complex<da::DAVector<T>>& c : l) same_env(e, c);
}

// Checks C++ only asserts (off in Release builds).
template <class T>
void check_base(const da::DAVector<T>& v, unsigned id) {
    if (id >= v.env_->layout().num_vars()) throw std::out_of_range("base id out of range");
}

// The engine resets the outputs before reading the inputs, so an output that
// is also an input gives a wrong result.
template <class T>
void check_disjoint(const std::vector<da::DAVector<T>>& out,
                    const std::vector<da::DAVector<T>>& in) {
    for (const da::DAVector<T>& o : out)
        for (const da::DAVector<T>& i : in)
            if (o.slot_ == i.slot_)
                throw std::invalid_argument("an output vector is also an input vector");
}

inline std::vector<int> exponents(const std::vector<int>& exps) {
    for (int e : exps)
        if (e < 0) throw std::invalid_argument("exponents must be non-negative");
    return exps;
}

inline void register_exceptions(nb::module_& m) {
    nb::exception<EnvError>(m, "EnvError", PyExc_RuntimeError);
    m.attr("EnvError").attr("__doc__") =
        "A DA vector was used after its env was closed or cleared, vectors of\n"
        "different envs met in one operation, or there is no current env.";
    nb::register_exception_translator([](const std::exception_ptr& p, void* payload) {
        try {
            std::rethrow_exception(p);
        } catch (const std::invalid_argument& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const std::domain_error& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const std::out_of_range& e) {
            PyErr_SetString(PyExc_IndexError, e.what());
        } catch (const std::logic_error& e) {
            PyErr_SetString(static_cast<PyObject*>(payload), e.what());
        }
    }, m.attr("EnvError").ptr());
}

void bind_env(nb::module_& m);
nb::class_<NDA> bind_nda(nb::module_& m);
void bind_cnda(nb::module_& m, nb::class_<NDA>& nda);
#ifdef DA_WITH_SYMBOLIC
void bind_expr(nb::module_& m);
nb::class_<SDA> bind_sda(nb::module_& m, nb::class_<NDA>& nda);
void bind_csda(nb::module_& m);
#endif
