// common.h — shared helpers of the C API (plan A.3, A.4).
#pragma once

#include "miradac.h"

#include "da/da.h"

#include <complex>
#include <stdexcept>
#include <vector>

namespace mdac {

using da::NDA;
using CNDA = std::complex<NDA>;
using NDAList = std::vector<NDA>;
using CNDAList = std::vector<CNDA>;

struct EnvError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Selects a DA object's env for the duration of a C API call.
struct EnvGuard {
    da::DAEnv* prev;
    explicit EnvGuard(da::DAEnv* e) {
        if (!e) throw EnvError("no DA environment");
        if (e->retired()) throw EnvError("DA environment has been cleared");
        prev = da::da_exchange_env(e);
    }
    ~EnvGuard() { da::da_exchange_env(prev); }
    EnvGuard(const EnvGuard&) = delete;
    EnvGuard& operator=(const EnvGuard&) = delete;
};

// Handles are the C++ objects themselves.
inline da::DAEnv* env(mdac_env* e) { return reinterpret_cast<da::DAEnv*>(e); }
inline mdac_env* handle(da::DAEnv* e) { return reinterpret_cast<mdac_env*>(e); }

// An env given as an argument, which must exist and not be retired.
inline da::DAEnv& live(mdac_env* h) {
    da::DAEnv* e = env(h);
    if (!e) throw EnvError("no DA environment");
    if (e->retired()) throw EnvError("DA environment has been cleared");
    return *e;
}
inline NDA& ref(mdac_nda* v) { return *reinterpret_cast<NDA*>(v); }
inline const NDA& ref(const mdac_nda* v) { return *reinterpret_cast<const NDA*>(v); }
inline mdac_nda* handle(NDA* v) { return reinterpret_cast<mdac_nda*>(v); }
inline NDAList& ref(mdac_ndalist* l) { return *reinterpret_cast<NDAList*>(l); }
inline const NDAList& ref(const mdac_ndalist* l) { return *reinterpret_cast<const NDAList*>(l); }
inline mdac_ndalist* handle(NDAList* l) { return reinterpret_cast<mdac_ndalist*>(l); }
inline CNDA& ref(mdac_cnda* v) { return *reinterpret_cast<CNDA*>(v); }
inline const CNDA& ref(const mdac_cnda* v) { return *reinterpret_cast<const CNDA*>(v); }
inline mdac_cnda* handle(CNDA* v) { return reinterpret_cast<mdac_cnda*>(v); }
inline CNDAList& ref(mdac_cndalist* l) { return *reinterpret_cast<CNDAList*>(l); }
inline const CNDAList& ref(const mdac_cndalist* l) { return *reinterpret_cast<const CNDAList*>(l); }
inline mdac_cndalist* handle(CNDAList* l) { return reinterpret_cast<mdac_cndalist*>(l); }

#ifdef DA_WITH_SYMBOLIC
using SymEngine::Expression;
using da::SDA;
using CSDA = std::complex<SDA>;
using SDAList = std::vector<SDA>;
using CSDAList = std::vector<CSDA>;

inline Expression& ref(mdac_expr* x) { return *reinterpret_cast<Expression*>(x); }
inline const Expression& ref(const mdac_expr* x) { return *reinterpret_cast<const Expression*>(x); }
inline mdac_expr* handle(Expression* x) { return reinterpret_cast<mdac_expr*>(x); }
inline SDA& ref(mdac_sda* v) { return *reinterpret_cast<SDA*>(v); }
inline const SDA& ref(const mdac_sda* v) { return *reinterpret_cast<const SDA*>(v); }
inline mdac_sda* handle(SDA* v) { return reinterpret_cast<mdac_sda*>(v); }
inline SDAList& ref(mdac_sdalist* l) { return *reinterpret_cast<SDAList*>(l); }
inline const SDAList& ref(const mdac_sdalist* l) { return *reinterpret_cast<const SDAList*>(l); }
inline mdac_sdalist* handle(SDAList* l) { return reinterpret_cast<mdac_sdalist*>(l); }
inline CSDA& ref(mdac_csda* v) { return *reinterpret_cast<CSDA*>(v); }
inline const CSDA& ref(const mdac_csda* v) { return *reinterpret_cast<const CSDA*>(v); }
inline mdac_csda* handle(CSDA* v) { return reinterpret_cast<mdac_csda*>(v); }
inline CSDAList& ref(mdac_csdalist* l) { return *reinterpret_cast<CSDAList*>(l); }
inline const CSDAList& ref(const mdac_csdalist* l) { return *reinterpret_cast<const CSDAList*>(l); }
inline mdac_csdalist* handle(CSDAList* l) { return reinterpret_cast<mdac_csdalist*>(l); }
#endif

// The engine indexes one pool with the slots of every vector it is given.
template <class T>
void check_same_env(const da::DAVector<T>& a, const da::DAVector<T>& b) {
    if (a.env_ != b.env_) throw EnvError("DA vectors belong to different environments");
}

// Exponents of a monomial; negative ones are invalid.
inline std::vector<int> exponents(const int* exps, size_t k) {
    for (size_t i = 0; i < k; ++i)
        if (exps[i] < 0) throw std::invalid_argument("exponents must be non-negative");
    return std::vector<int>(exps, exps + k);
}

// Maps the exception being handled to a status and records its message.
mdac_status current_exception_status() noexcept;

// MDAC_ERR_UNSUPPORTED, for a symbolic function in a numeric-only build.
mdac_status unsupported() noexcept;

} // namespace mdac

// No exception crosses the C boundary.
#define MDAC_TRY try
#define MDAC_CATCH                                                  \
    catch (...) { return mdac::current_exception_status(); }       \
    return MDAC_OK;

// The body of a symbolic function: the code itself, or, in a numeric-only
// build, `fallback` (MDAC_SYM: return MDAC_ERR_UNSUPPORTED).
#ifdef DA_WITH_SYMBOLIC
#define MDAC_SYM_ELSE(fallback, ...) __VA_ARGS__
#else
#define MDAC_SYM_ELSE(fallback, ...) fallback
#endif
#define MDAC_SYM(...) MDAC_SYM_ELSE(return mdac::unsupported();, __VA_ARGS__)
