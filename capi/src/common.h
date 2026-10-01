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

// The engine indexes one pool with the slots of every vector it is given.
template <class T>
void check_same_env(const da::DAVector<T>& a, const da::DAVector<T>& b) {
    if (a.env_ != b.env_) throw EnvError("DA vectors belong to different environments");
}

// Maps the exception being handled to a status and records its message.
mdac_status current_exception_status() noexcept;

} // namespace mdac

// No exception crosses the C boundary.
#define MDAC_TRY try
#define MDAC_CATCH                                                  \
    catch (...) { return mdac::current_exception_status(); }       \
    return MDAC_OK;
