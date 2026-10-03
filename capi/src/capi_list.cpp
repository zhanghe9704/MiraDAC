// capi_list.cpp — NDA lists and the NDA algorithms of the C API (plan T3.1).
#include "common.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <set>
#include <vector>

using namespace mdac;

namespace {

// The algorithms take non-const lists but do not change their inputs.
NDAList& mut(const NDAList& l) { return const_cast<NDAList&>(l); }

// Every vector of a list belongs to one env, so a list is freed in that env.
void check_joins(const NDAList& l, const NDA& v) {
    if (!l.empty()) check_same_env(l[0], v);
}

void same_env(const NDA& first, const NDAList& l) {
    for (const NDA& v : l) check_same_env(first, v);
}

// Checks C++ only asserts (python/src/common.h check_base).
void check_base(const NDA& v, unsigned id) {
    if (id >= v.env_->layout().num_vars()) throw std::out_of_range("base id out of range");
}

std::vector<unsigned> check_ids(const NDA& first, const unsigned* ids, size_t k, const NDAList& v) {
    if (k != v.size()) throw std::invalid_argument("base ids and vectors must have the same length");
    std::vector<unsigned> b(ids, ids + k);
    if (std::set<unsigned>(b.begin(), b.end()).size() != b.size())
        throw std::invalid_argument("duplicate base id");
    for (unsigned id : b) check_base(first, id);
    same_env(first, v);
    return b;
}

void check_point(const NDAList& m, size_t n) {
    if (n != m[0].env_->layout().num_vars())
        throw std::invalid_argument("the point must have nvars coordinates");
    same_env(m[0], m);
}

// A new list of n zero vectors in the current env.
mdac_ndalist* zeros(size_t n) {
    auto* l = new NDAList();
    try {
        l->reserve(n);
        for (size_t i = 0; i < n; ++i) l->emplace_back();
    } catch (...) {
        delete l;
        throw;
    }
    return handle(l);
}

// The kernel of the Python evaluate_map (python/src/bind_nda.cpp eval_points),
// copied: out[p * k + c] = map[c] at point p, for k = coef.size(). Points go
// in blocks of kB with the point index innermost, so the compiler vectorizes
// over points; each output still sums its terms in ascending monomial order,
// as ad_composition does, so the results equal compose(map, point) bit for bit.
// GCC's target_clones IFUNC resolver is unsupported by musl.
#if defined(__GNUC__) && defined(__x86_64__) && defined(__linux__) && defined(__GLIBC__)
__attribute__((target_clones("avx2", "default")))
#endif
void eval_points(const double* pts, size_t n, size_t nv, size_t nd,
                 const std::vector<unsigned>& exps, const std::vector<const double*>& coef,
                 const std::vector<unsigned>& lens, double* buf) {
    const size_t k = coef.size(), len = exps.size() / nv;
    constexpr size_t kB = 128;
    std::vector<double> power(nv * (nd + 1) * kB), acc(k * kB);
    double prod[kB];
    for (size_t b = 0; b < n; b += kB) {
        const size_t nb = std::min(kB, n - b);
        const double* x = pts + b * nv;
        for (size_t j = 0; j < nv; ++j) {
            double* pw = power.data() + j * (nd + 1) * kB;
            for (size_t q = 0; q < nb; ++q) pw[q] = 1.0;
            for (size_t o = 1; o <= nd; ++o)
                for (size_t q = 0; q < nb; ++q)
                    pw[o * kB + q] = x[q * nv + j] * pw[(o - 1) * kB + q];
        }
        for (size_t c = 0; c < k; ++c)
            for (size_t q = 0; q < nb; ++q) acc[c * kB + q] = coef[c][0];
        for (size_t i = 1; i < len; ++i) {
            const unsigned* ex = exps.data() + i * nv;
            std::memcpy(prod, power.data() + ex[0] * kB, nb * sizeof(double));
            for (size_t j = 1; j < nv; ++j) {
                const double* __restrict pw = power.data() + (j * (nd + 1) + ex[j]) * kB;
                for (size_t q = 0; q < nb; ++q) prod[q] *= pw[q];
            }
            for (size_t c = 0; c < k; ++c) {
                if (i >= lens[c]) continue;
                const double a = coef[c][i];
                double* __restrict ac = acc.data() + c * kB;
                for (size_t q = 0; q < nb; ++q) ac[q] += prod[q] * a;
            }
        }
        for (size_t q = 0; q < nb; ++q)
            for (size_t c = 0; c < k; ++c) buf[(b + q) * k + c] = acc[c * kB + q];
    }
}

} // namespace

extern "C" {

// ---- Lists ------------------------------------------------------------------

mdac_status mdac_ndalist_new(mdac_ndalist** out) {
    MDAC_TRY { *out = handle(new NDAList()); } MDAC_CATCH
}

mdac_status mdac_ndalist_from(const mdac_nda* const* vs, size_t n, mdac_ndalist** out) {
    MDAC_TRY {
        if (n == 0) { *out = handle(new NDAList()); return MDAC_OK; }
        EnvGuard g(ref(vs[0]).env_);
        for (size_t i = 1; i < n; ++i) check_same_env(ref(vs[0]), ref(vs[i]));
        auto* l = new NDAList();
        try {
            l->reserve(n);
            for (size_t i = 0; i < n; ++i) l->push_back(ref(vs[i]));
        } catch (...) {
            delete l;
            throw;
        }
        *out = handle(l);
    } MDAC_CATCH
}

void mdac_ndalist_free(mdac_ndalist* l) { delete &ref(l); }

size_t mdac_ndalist_length(const mdac_ndalist* l) { return ref(l).size(); }

mdac_env* mdac_ndalist_env(const mdac_ndalist* l) {
    return ref(l).empty() ? nullptr : handle(ref(l)[0].env_);
}

mdac_status mdac_ndalist_get(const mdac_ndalist* l_, size_t i, mdac_nda** out) {
    MDAC_TRY {
        const NDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(l[i].env_);
        *out = handle(new NDA(l[i]));
    } MDAC_CATCH
}

mdac_status mdac_ndalist_set(mdac_ndalist* l_, size_t i, const mdac_nda* v) {
    MDAC_TRY {
        NDAList& l = ref(l_);
        if (i >= l.size()) throw std::out_of_range("list index out of range");
        EnvGuard g(ref(v).env_);
        check_same_env(l[i], ref(v));
        l[i] = ref(v);
    } MDAC_CATCH
}

mdac_status mdac_ndalist_push(mdac_ndalist* l_, const mdac_nda* v) {
    MDAC_TRY {
        NDAList& l = ref(l_);
        EnvGuard g(ref(v).env_);
        check_joins(l, ref(v));
        l.push_back(ref(v));
    } MDAC_CATCH
}

// ---- Algorithms -------------------------------------------------------------

mdac_status mdac_nda_der(const mdac_nda* v, unsigned i, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        check_base(ref(v), i);
        *out = handle(new NDA(da::da_der(ref(v), i)));
    } MDAC_CATCH
}

mdac_status mdac_nda_integ(const mdac_nda* v, unsigned i, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(v).env_);
        check_base(ref(v), i);
        *out = handle(new NDA(da::da_int(ref(v), i)));
    } MDAC_CATCH
}

mdac_status mdac_nda_substitute_d(const mdac_nda* iv, unsigned i, double x, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(iv).env_);
        check_base(ref(iv), i);
        NDA o;
        da::da_substitute_const(ref(iv), i, x, o);
        *out = handle(new NDA(std::move(o)));
    } MDAC_CATCH
}

mdac_status mdac_nda_substitute(const mdac_nda* iv, unsigned i, const mdac_nda* v, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(iv).env_);
        check_same_env(ref(iv), ref(v));
        check_base(ref(iv), i);
        NDA o;
        da::da_substitute(ref(iv), i, ref(v), o);
        *out = handle(new NDA(std::move(o)));
    } MDAC_CATCH
}

mdac_status mdac_nda_substitute_multi(const mdac_nda* iv, const unsigned* ids, size_t k,
                                      const mdac_ndalist* v, mdac_nda** out) {
    MDAC_TRY {
        EnvGuard g(ref(iv).env_);
        std::vector<unsigned> b = check_ids(ref(iv), ids, k, ref(v));
        NDA o;
        da::da_substitute(ref(iv), b, mut(ref(v)), o);
        *out = handle(new NDA(std::move(o)));
    } MDAC_CATCH
}

mdac_status mdac_ndalist_substitute(const mdac_ndalist* m_, const unsigned* ids, size_t k,
                                    const mdac_ndalist* v, mdac_ndalist** out) {
    MDAC_TRY {
        const NDAList& m = ref(m_);
        if (m.empty()) { *out = handle(new NDAList()); return MDAC_OK; }
        EnvGuard g(m[0].env_);
        std::vector<unsigned> b = check_ids(m[0], ids, k, ref(v));
        same_env(m[0], m);
        mdac_ndalist* o = zeros(m.size());
        try {
            da::da_substitute(mut(m), b, mut(ref(v)), ref(o));
        } catch (...) {
            mdac_ndalist_free(o);
            throw;
        }
        *out = o;
    } MDAC_CATCH
}

mdac_status mdac_ndalist_compose(const mdac_ndalist* m_, const mdac_ndalist* v_, mdac_ndalist** out) {
    MDAC_TRY {
        const NDAList &m = ref(m_), &v = ref(v_);
        if (m.empty()) { *out = handle(new NDAList()); return MDAC_OK; }
        EnvGuard g(m[0].env_);
        if (v.size() != m[0].env_->layout().num_vars())
            throw std::invalid_argument("compose: the arguments must be nvars vectors");
        same_env(m[0], m);
        same_env(m[0], v);
        mdac_ndalist* o = zeros(m.size());
        try {
            da::da_composition(mut(m), mut(v), ref(o));
        } catch (...) {
            mdac_ndalist_free(o);
            throw;
        }
        *out = o;
    } MDAC_CATCH
}

mdac_status mdac_ndalist_compose_d(const mdac_ndalist* m_, const double* pt, size_t n, double* out) {
    MDAC_TRY {
        const NDAList& m = ref(m_);
        if (m.empty()) return MDAC_OK;
        EnvGuard g(m[0].env_);
        check_point(m, n);
        std::vector<double> p(pt, pt + n), o(m.size());
        da::da_composition(mut(m), p, o);
        std::copy(o.begin(), o.end(), out);
    } MDAC_CATCH
}

mdac_status mdac_ndalist_compose_z(const mdac_ndalist* m_, const double* pt, size_t n, double* out) {
    MDAC_TRY {
        const NDAList& m = ref(m_);
        if (m.empty()) return MDAC_OK;
        EnvGuard g(m[0].env_);
        check_point(m, n);
        const auto* z = reinterpret_cast<const std::complex<double>*>(pt);
        std::vector<std::complex<double>> p(z, z + n), o(m.size());
        da::da_composition(mut(m), p, o);
        std::copy(o.begin(), o.end(), reinterpret_cast<std::complex<double>*>(out));
    } MDAC_CATCH
}

mdac_status mdac_ndalist_inv_map(const mdac_ndalist* m_, int dim, mdac_ndalist** out) {
    MDAC_TRY {
        const NDAList& m = ref(m_);
        if (m.empty()) throw std::invalid_argument("inv_map: empty map");
        EnvGuard g(m[0].env_);
        if (dim < 1 || dim > static_cast<int>(m[0].env_->layout().num_vars()))
            throw std::invalid_argument("inv_map: dim must be in [1, nvars]");
        if (m.size() < static_cast<size_t>(dim))
            throw std::invalid_argument("inv_map: the map needs at least dim vectors");
        same_env(m[0], m);
        for (const NDA& v : m)
            if (std::fabs(v.con()) >= std::numeric_limits<double>::min())
                throw std::invalid_argument("inv_map: the map must have zero constant parts");
        mdac_ndalist* o = zeros(static_cast<size_t>(dim));
        try {
            da::inv_map(mut(m), dim, ref(o));
        } catch (...) {
            mdac_ndalist_free(o);
            throw;
        }
        *out = o;
    } MDAC_CATCH
}

mdac_status mdac_ndalist_evaluate_map(const mdac_ndalist* m_, const double* pts, size_t npts,
                                      double* out) {
    MDAC_TRY {
        const NDAList& m = ref(m_);
        if (m.empty()) return MDAC_OK;
        da::DAEnv* e = m[0].env_;
        EnvGuard g(e);
        same_env(m[0], m);
        da::Layout& l = e->layout();
        da::Pool<double>& pool = e->pool<double>();
        const size_t nv = l.num_vars(), nd = l.max_order();

        size_t len = 0;
        std::vector<unsigned> lens;
        std::vector<const double*> coef;
        for (const NDA& v : m) {
            lens.push_back(pool.len(v.slot_));
            coef.push_back(pool.slot(v.slot_));
            len = std::max<size_t>(len, lens.back());
        }
        std::vector<unsigned> exps(len * nv);
        for (size_t i = 1; i < len; ++i) {
            const std::vector<int>& o = l.orders_ref(static_cast<unsigned>(i));
            std::copy(o.begin(), o.begin() + nv, exps.begin() + i * nv);
        }
        eval_points(pts, npts, nv, nd, exps, coef, lens, out);
    } MDAC_CATCH
}

mdac_status mdac_nda_eval(const mdac_nda* v_, const double* pt, size_t n, double* out) {
    MDAC_TRY {
        const NDA& v = ref(v_);
        EnvGuard g(v.env_);
        da::Layout& l = v.env_->layout();
        if (n != l.num_vars()) throw std::invalid_argument("the point must have nvars coordinates");
        std::vector<unsigned> iv{v.slot_};
        std::vector<double> p(pt, pt + n), o(1);
        da::detail::ad_composition(l, v.env_->pool<double>(), iv, p, o);
        *out = o[0];
    } MDAC_CATCH
}

mdac_status mdac_exponents(mdac_env* e_, int* buf, size_t cap, size_t* n) {
    MDAC_TRY {
        da::DAEnv* e = env(e_);
        EnvGuard g(e);
        da::Layout& l = e->layout();
        const size_t len = l.full_len(), nv = l.num_vars();
        *n = len * nv;
        for (size_t i = 0; i < len && (i + 1) * nv <= cap; ++i) {
            const std::vector<int>& o = l.orders_ref(static_cast<unsigned>(i));
            std::copy(o.begin(), o.begin() + nv, buf + i * nv);
        }
    } MDAC_CATCH
}

}
