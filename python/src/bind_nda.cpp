// bind_nda.cpp — NDA class, NDAList, var(), the NDA math functions and the
// NDA algorithms (plan T1.3-T1.5, T2.1-T2.4).
#include "arith.h"

#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <sstream>

namespace {

using Points = nb::ndarray<const double, nb::ndim<2>, nb::c_contig>;
using Complexes = std::vector<std::complex<double>>;

// Returns a heap buffer as a numpy array that owns it.
template <class T, class... Shape>
nb::ndarray<nb::numpy, T> owned(T* buf, Shape... shape) {
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<T*>(p); });
    return nb::ndarray<nb::numpy, T>(buf, {static_cast<size_t>(shape)...}, owner);
}

std::vector<unsigned> slots(const NDAList& l) {
    std::vector<unsigned> s;
    s.reserve(l.size());
    for (const NDA& v : l) s.push_back(v.slot_);
    return s;
}

template <NDA (*F)(const NDA&)>
auto unary(const NDA& v) {
    EnvGuard g(v.env_);
    return make_da<NDA>([&] { return F(v); });
}

// evaluate_map's kernel: out[p * k + c] = map[c] at point p, for k = coef.size().
// Points go in blocks of kB, with the point index innermost so the compiler
// vectorizes over points; each output still sums its terms in ascending
// monomial order, as ad_composition does. AVX2 without FMA contraction
// rounds exactly like the default build.
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
            // 1.0 * pw[q] == pw[q], so starting from the first factor
            // keeps the rounding of ad_composition's product.
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

// T2.1-T2.4.
static void bind_algorithms(nb::module_& m) {
    // nb::bind_vector registers implicitly_convertible<iterable, NDAList>, so
    // a plain Python list is accepted wherever an NDAList& is expected; it is
    // copied into a temporary NDAList, one NDA copy per element. A temporary
    // cannot carry results back, so list outputs take .noconvert().
    nb::bind_vector<NDAList>(m, "NDAList",
        "A list of NDA held in C++, for the map-level functions (no per-element\n"
        "copies). A plain list is accepted as an input too, copied element by\n"
        "element; output arguments must be an NDAList.");

    m.def("da_der", [](const NDA& v, unsigned id) {
        EnvGuard g(v.env_);
        check_base(v, id);
        return da::da_der(v, id);
    }, nb::arg("v"), nb::arg("base_id"),
       "Derivative of v with respect to variable base_id (alias der). With\n"
       "result given, it is written there and None is returned.");
    m.def("da_der", [](const NDA& v, unsigned id, NDA& result) {
        EnvGuard g(v.env_);
        same_env(v.env_, result);
        check_base(v, id);
        da::da_der(v, id, result);
    }, nb::arg("v"), nb::arg("base_id"), nb::arg("result"));
    m.def("da_int", [](const NDA& v, unsigned id) {
        EnvGuard g(v.env_);
        check_base(v, id);
        return da::da_int(v, id);
    }, nb::arg("v"), nb::arg("base_id"),
       "Integral of v with respect to variable base_id (alias int_). With\n"
       "result given, it is written there and None is returned.");
    m.def("da_int", [](const NDA& v, unsigned id, NDA& result) {
        EnvGuard g(v.env_);
        same_env(v.env_, result);
        check_base(v, id);
        da::da_int(v, id, result);
    }, nb::arg("v"), nb::arg("base_id"), nb::arg("result"));

    m.def("da_substitute_const", [](const NDA& iv, unsigned id, double x, NDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute_const(iv, id, x, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("x"), nb::arg("ov"),
       "ov = iv with the number x substituted for variable base_id.");
    m.def("da_substitute", [](const NDA& iv, unsigned id, const NDA& v, NDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, v);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute(iv, id, v, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("v"), nb::arg("ov"),
       "Substitute into iv and write the result to ov (alias substitute).\n\n"
       "base_id is one variable and v a DA vector or a number, or base_id is a\n"
       "list of variables and v a list of DA vectors, one per variable. iv and\n"
       "ov may also be lists (ov then an NDAList or SDAList of len(iv)).");
    m.def("da_substitute", [](const NDA& iv, unsigned id, double x, NDA& ov) {
        EnvGuard g(iv.env_);
        same_env(iv.env_, ov);
        check_base(iv, id);
        da::da_substitute(iv, id, x, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("x"), nb::arg("ov"));

    // Checks shared by the multi-base forms; C++ only asserts them.
    auto check_multi = [](const NDA& first, const std::vector<unsigned>& ids,
                          const NDAList& v) {
        if (ids.size() != v.size())
            throw std::invalid_argument("base_id and v must have the same length");
        if (std::set<unsigned>(ids.begin(), ids.end()).size() != ids.size())
            throw std::invalid_argument("duplicate base id");
        for (unsigned id : ids) check_base(first, id);
        same_env(first.env_, v);
    };
    m.def("da_substitute", [check_multi](const NDA& iv, std::vector<unsigned> ids,
                                         NDAList& v, NDA& ov) {
        EnvGuard g(iv.env_);
        check_multi(iv, ids, v);
        same_env(iv.env_, ov);
        da::da_substitute(iv, ids, v, ov);
    }, nb::arg("iv"), nb::arg("base_id"), nb::arg("v"), nb::arg("ov"));
    m.def("da_substitute", [check_multi](NDAList& ivecs, std::vector<unsigned> ids,
                                         NDAList& v, NDAList& ovecs) {
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

    // The float and complex forms return the values: a Python list cannot be
    // an output argument.
    auto check_point = [](const NDAList& ivecs, size_t n) {
        if (n != ivecs[0].env_->layout().num_vars())
            throw std::invalid_argument("v must have nvars() entries");
        same_env(ivecs[0].env_, ivecs);
    };
    m.def("da_composition", [](NDAList& ivecs, NDAList& v, NDAList& ovecs) {
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
    }, nb::arg("ivecs"), nb::arg("v"), nb::arg("ovecs").noconvert(),
       "Compose a map with v (alias compose): ivecs[i](v[0], ..., v[nvars-1]).\n\n"
       "With DA vectors v, the results go to ovecs (a list of len(ivecs)) and\n"
       "None is returned; with numbers (float or complex) v, the list of\n"
       "values is returned. v has nvars() entries.");
    m.def("da_composition", [check_point](NDAList& ivecs, std::vector<double> v) {
        std::vector<double> out(ivecs.size());
        if (ivecs.empty()) return out;
        EnvGuard g(ivecs[0].env_);
        check_point(ivecs, v.size());
        da::da_composition(ivecs, v, out);
        return out;
    }, nb::arg("ivecs"), nb::arg("v"));
    m.def("da_composition", [check_point](NDAList& ivecs, Complexes v) {
        Complexes out(ivecs.size());
        if (ivecs.empty()) return out;
        EnvGuard g(ivecs[0].env_);
        check_point(ivecs, v.size());
        da::da_composition(ivecs, v, out);
        return out;
    }, nb::arg("ivecs"), nb::arg("v"));

    m.def("inv_map", [](NDAList& ivecs, int dim, NDAList& ovecs) {
        if (ivecs.empty()) throw std::invalid_argument("inv_map: empty map");
        EnvGuard g(ivecs[0].env_);
        if (dim < 1 || dim > static_cast<int>(ivecs[0].env_->layout().num_vars()))
            throw std::invalid_argument("inv_map: dim must be in [1, nvars()]");
        if (ivecs.size() < static_cast<size_t>(dim) || ovecs.size() < static_cast<size_t>(dim))
            throw std::invalid_argument("inv_map: ivecs and ovecs need at least dim vectors");
        same_env(ivecs[0].env_, ivecs);
        same_env(ivecs[0].env_, ovecs);
        for (const NDA& v : ivecs)
            if (std::fabs(v.con()) >= std::numeric_limits<double>::min())
                throw std::invalid_argument("inv_map: the map must have zero constant parts");
        da::inv_map(ivecs, dim, ovecs);
    }, nb::arg("ivecs"), nb::arg("dim"), nb::arg("ovecs").noconvert(),
       "Inverse of the map given by the first dim vectors of ivecs, written to\n"
       "the first dim vectors of ovecs (an NDAList). The map must have zero\n"
       "constant parts.");

    m.def("devide_by_element", [](NDA& t, NDA& b) {
        EnvGuard g(t.env_);
        same_env(t.env_, b);
        return da::devide_by_element(t, b);
    }, nb::arg("t"), nb::arg("b"),
       "Coefficient-wise t / b (a coefficient of t stays where b's is zero).");
    m.def("read_da_from_file", [](const std::string& filename, NDA& d) {
        EnvGuard g(d.env_);
        return da::read_da_from_file(filename, d);
    }, nb::arg("filename"), nb::arg("d"),
       "Read an NDA written by str()/print into d; False if the file is\n"
       "incomplete.");
    m.def("compare_da_vectors", [](NDA& a, NDA& b, double eps) {
        EnvGuard g(a.env_);
        return da::compare_da_vectors(a, b, eps);
    }, nb::arg("a"), nb::arg("b"), nb::arg("eps") = 1e-15,
       "True if a == b to eps: every coefficient of a - b is below eps times\n"
       "b's (below eps where b's is zero).");

    m.attr("der") = m.attr("da_der");
    m.attr("int_") = m.attr("da_int");
    m.attr("substitute") = m.attr("da_substitute");
    m.attr("compose") = m.attr("da_composition");

    // Evaluates like the double form of da::detail::ad_composition (same
    // products, same summation order, so at finite points the results equal
    // da_composition(map, point) bit for bit), but does not call it: a loop
    // calling it, which rebuilds its power tables on every call, was only
    // ~7x faster than a Python loop, and T2.3 asks for 20x.
    m.def("evaluate_map", [](NDAList& map, Points pts) {
        const size_t n = pts.shape(0), k = map.size();
        double* buf = new double[n * k];
        auto result = owned(buf, n, k);
        if (k == 0) return result;
        da::DAEnv* e = map[0].env_;
        EnvGuard g(e);
        same_env(e, map);
        da::Layout& l = e->layout();
        da::Pool<double>& pool = e->pool<double>();
        const size_t nv = l.num_vars(), nd = l.max_order();
        if (pts.shape(1) != nv) throw std::invalid_argument("pts must have shape (N, nvars())");

        size_t len = 0;
        std::vector<unsigned> lens;
        std::vector<const double*> coef;
        for (const NDA& v : map) {
            lens.push_back(pool.len(v.slot_));
            coef.push_back(pool.slot(v.slot_));
            len = std::max<size_t>(len, lens.back());
        }
        std::vector<unsigned> exps(len * nv);
        for (size_t i = 1; i < len; ++i) {
            const std::vector<int>& o = l.orders_ref(static_cast<unsigned>(i));
            std::copy(o.begin(), o.begin() + nv, exps.begin() + i * nv);
        }

        eval_points(pts.data(), n, nv, nd, exps, coef, lens, buf);
        return result;
    }, nb::arg("map"), nb::arg("pts"),
       "Evaluate a map at many points in one C++ loop.\n\n"
       "pts is a float64 array of shape (N, nvars()); the result has shape\n"
       "(N, len(map)) and equals compose(map, pts[p]) for each point p.");

    m.def("exponents", [](da::DAEnv* env) {
        da::DAEnv* e = env_or_current(env);
        EnvGuard g(e);
        da::Layout& l = e->layout();
        const size_t n = l.full_len(), nv = l.num_vars();
        int32_t* buf = new int32_t[n * nv];
        for (size_t i = 0; i < n; ++i) {
            const std::vector<int>& o = l.orders_ref(static_cast<unsigned>(i));
            std::copy(o.begin(), o.begin() + nv, buf + i * nv);
        }
        return owned(buf, n, nv);
    }, nb::arg("env") = nb::none(),
       "The exponents of every monomial, an int32 array of shape\n"
       "(full_length(), nvars()): row i belongs to coefficient i. For the\n"
       "current env, or env.");
}

nb::class_<NDA> bind_nda(nb::module_& m) {
    nb::class_<NDA> c(m, "NDA",
        "A numeric DA vector: a truncated power series in nvars variables with\n"
        "float coefficients. It lives in the env it was created in.");
    c.def("__init__", [](NDA* self, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) NDA();
    }, nb::kw_only(), nb::arg("env") = nb::none(),
       "NDA(): zero. NDA(x): the constant x. In the current env, or env.");
    c.def("__init__", [](NDA* self, double x, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        new (self) NDA(x);
    }, nb::arg("x"), nb::kw_only(), nb::arg("env") = nb::none());
    c.def_static("from_coeffs", [](nb::ndarray<const double, nb::ndim<1>, nb::c_contig> a,
                                   da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        NDA v;
        da::Pool<double>& pool = v.env_->pool<double>();
        if (a.shape(0) > v.env_->layout().full_len())
            throw std::invalid_argument("from_coeffs: more coefficients than full_length()");
        const unsigned n = static_cast<unsigned>(a.shape(0));
        std::memcpy(pool.slot(v.slot_), a.data(), n * sizeof(double));
        pool.set_len(v.slot_, n);
        return v;
    }, nb::arg("a"), nb::kw_only(), nb::arg("env") = nb::none(),
       "An NDA with the coefficients a (float64, at most full_length()\n"
       "values, in monomial order; see exponents()).");

    c.def("copy", [](const NDA& v) { EnvGuard g(v.env_); return NDA(v); }, "A copy.");
    c.def("__copy__", [](const NDA& v) { EnvGuard g(v.env_); return NDA(v); }, "A copy.");
    c.def("__deepcopy__", [](const NDA& v, nb::handle) { EnvGuard g(v.env_); return NDA(v); },
          nb::arg("memo"), "A copy.");

    c.def_prop_rw("con",
        [](const NDA& v) { EnvGuard g(v.env_); return v.con(); },
        [](NDA& v, double x) { EnvGuard g(v.env_); v.reset_const(x); },
        "The constant part. Setting it resets the vector to that constant.");
    c.def_prop_ro("length", [](const NDA& v) { EnvGuard g(v.env_); return v.length(); },
                  "Number of stored coefficients (up to the last non-zero one).");
    c.def_prop_ro("n_element", [](const NDA& v) { EnvGuard g(v.env_); return v.n_element(); },
                  "Number of non-zero coefficients.");
    c.def_prop_ro("nvars", [](const NDA& v) {
        EnvGuard g(v.env_);
        return v.env_->layout().num_vars();
    }, "Number of variables of the vector's env.");
    c.def_prop_ro("order", [](const NDA& v) {
        EnvGuard g(v.env_);
        return v.env_->layout().max_order();
    }, "Current truncation order of the vector's env.");
    c.def("norm", [](const NDA& v) { EnvGuard g(v.env_); return v.norm(); },
          "The largest coefficient magnitude.");
    c.def("weighted_norm", [](const NDA& v, double w) {
        EnvGuard g(v.env_);
        return v.weighted_norm(w);
    }, nb::arg("w"), "The largest magnitude of a coefficient times w**(its order).");
    c.def("iszero", [](const NDA& v, std::optional<double> eps) {
        EnvGuard g(v.env_);
        return eps ? v.iszero(*eps) : v.iszero();
    }, nb::arg("eps") = nb::none(),
       "True if every coefficient is below eps in magnitude (default get_eps()).");
    c.def("clean", [](NDA& v, std::optional<double> eps) {
        EnvGuard g(v.env_);
        eps ? v.clean(*eps) : v.clean();
    }, nb::arg("eps") = nb::none(),
       "Set coefficients below eps in magnitude (default get_eps()) to zero.");
    c.def("reset", [](NDA& v) { EnvGuard g(v.env_); v.reset(); }, "Set every coefficient to zero.");

    c.def("element", [](NDA& v, int i) { EnvGuard g(v.env_); return v.element(i); },
          nb::arg("i"),
          "element(i): coefficient i in monomial order (see exponents()), 0 past\n"
          "length. element(exps): the coefficient of the monomial with the given\n"
          "exponent of each variable.");
    c.def("element", [](NDA& v, const std::vector<int>& exps) {
        EnvGuard g(v.env_);
        return v.element(exponents(exps));
    }, nb::arg("exps"));
    c.def("set_element", [](NDA& v, const std::vector<int>& exps, double value) {
        EnvGuard g(v.env_);
        v.set_element(exponents(exps), value);
    }, nb::arg("exps"), nb::arg("value"),
       "Set the coefficient of the monomial with exponents exps.");
    c.def("index_element", [](const NDA& v, unsigned i) {
        EnvGuard g(v.env_);
        std::vector<unsigned> exps;
        double value;
        v.element(i, exps, value);
        return std::make_pair(std::vector<int>(exps.begin(), exps.end()), value);
    }, nb::arg("i"), "(exponents, coefficient) of coefficient i in monomial order.");
    c.def("coeffs", [](const NDA& v) {
        EnvGuard g(v.env_);
        const size_t n = v.length();
        double* buf = new double[n];
        std::memcpy(buf, v.env_->pool<double>().slot(v.slot_), n * sizeof(double));
        return owned(buf, n);
    }, "A copy of the length stored coefficients, a float64 array.");
    c.def("__call__", [](const NDA& v, const std::vector<double>& pt) {
        EnvGuard g(v.env_);
        da::Layout& l = v.env_->layout();
        if (pt.size() != l.num_vars())
            throw std::invalid_argument("point must have nvars() coordinates");
        std::vector<unsigned> iv{v.slot_};
        std::vector<double> p(pt), out(1);
        da::detail::ad_composition(l, v.env_->pool<double>(), iv, p, out);
        return out[0];
    }, nb::arg("pt"), "The value at the point pt (nvars coordinates).");

    c.def("__repr__", [](const NDA& v) {
        EnvGuard g(v.env_);
        const da::Layout& l = v.env_->layout();
        std::ostringstream os;
        os << "NDA(order=" << l.max_order() << ", nvars=" << l.num_vars()
           << ", nonzero=" << v.n_element() << ")";
        return os.str();
    }, "NDA(order=..., nvars=..., nonzero=...).");
    c.def("__str__", [](const NDA& v) {
        EnvGuard g(v.env_);
        std::ostringstream os;
        os << v;
        return os.str();
    }, "The table of terms that C++ prints.");

    bind_arith(c);
    bind_env_property(c);

    m.def("var", [](int i, da::DAEnv* env) {
        EnvGuard g(env_or_current(env));
        if (i < 0 || i >= static_cast<int>(da::da_current_env().layout().num_vars()))
            throw std::out_of_range("var: index out of range");
        return da::da_base(static_cast<unsigned>(i));
    }, nb::arg("i"), nb::kw_only(), nb::arg("env") = nb::none(),
       "The base vector of variable i (0 <= i < nvars), in the current env or\n"
       "env. base[i] is the same.");

    m.def("sqrt", &unary<da::sqrt>, "Square root (NDA, CNDA, SDA or CSDA).");
    m.def("exp", &unary<da::exp>, "Exponential (NDA, CNDA, SDA or CSDA).");
    m.def("log", &unary<da::log>, "Natural logarithm (NDA, CNDA, SDA or CSDA).");
    m.def("sin", &unary<da::sin>, "Sine (NDA or SDA).");
    m.def("cos", &unary<da::cos>, "Cosine (NDA or SDA).");
    m.def("tan", &unary<da::tan>, "Tangent (NDA or SDA).");
    m.def("asin", &unary<da::asin>, "Arcsine (NDA, CNDA, SDA or CSDA).");
    m.def("acos", &unary<da::acos>, "Arccosine (NDA, CNDA, SDA or CSDA).");
    m.def("atan", &unary<da::atan>, "Arctangent (NDA, CNDA, SDA or CSDA).");
    m.def("sinh", &unary<da::sinh>, "Hyperbolic sine (NDA or SDA).");
    m.def("cosh", &unary<da::cosh>, "Hyperbolic cosine (NDA or SDA).");
    m.def("tanh", &unary<da::tanh>, "Hyperbolic tangent (NDA or SDA).");
    m.def("asinh", &unary<da::asinh>, "Inverse hyperbolic sine (NDA, CNDA or CSDA).");
    m.def("acosh", &unary<da::acosh>, "Inverse hyperbolic cosine (NDA, CNDA or CSDA).");
    m.def("atanh", &unary<da::atanh>, "Inverse hyperbolic tangent (NDA, CNDA or CSDA).");
    m.def("erf", &unary<da::erf>, "Error function (NDA or SDA).");
    m.def("abs", [](const NDA& v) {
        EnvGuard g(v.env_);
        return da::abs(v);
    }, "abs(NDA): the norm(), a float.");
    m.def("pow", [](const NDA& v, int n) { EnvGuard g(v.env_); return da::pow(v, n); },
          "v raised to an int or float power (NDA, CNDA, SDA or CSDA); the same\n"
          "as v ** p.");
    m.def("pow", [](const NDA& v, double x) { EnvGuard g(v.env_); return da::pow(v, x); });

    m.def("compare_da_with_file", [](const std::string& filename, NDA& d, double eps) {
        EnvGuard g(d.env_);
        return da::compare_da_with_file(filename, d, eps);
    }, nb::arg("filename"), nb::arg("d"), nb::arg("eps") = 1e-15,
       "compare_da_vectors of the NDA stored in filename and d.");

    bind_algorithms(m);
    return c;
}
