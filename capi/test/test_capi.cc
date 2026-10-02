// test_capi.cc — the C API against direct C++ results (plan T1.2-T7.1).
//
// The C API library holds its own hidden copy of the da library, so the C++
// reference vectors below live in a separate default env of the same shape.
#include "catch.hpp"
#include "miradac.h"

#include "da/da.h"

#include <cmath>
#include <complex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

#define OK(call) REQUIRE((call) == MDAC_OK)

constexpr unsigned kOrder = 4, kNvars = 3, kPool = 200;

// A C API default env and a C++ default env of the same shape.
struct Envs {
    mdac_env* env = nullptr;
    explicit Envs(unsigned pool = kPool) {
        OK(mdac_init(kOrder, kNvars, pool, 0));
        OK(mdac_env_default(&env));
        da::da_init(kOrder, kNvars, pool);
    }
    ~Envs() {
        mdac_clear();
        da::da_clear();
    }
};

// Owns a C API vector.
struct H {
    mdac_nda* p = nullptr;
    H() = default;
    explicit H(mdac_nda* q) : p(q) {}
    H(const H&) = delete;
    H& operator=(const H&) = delete;
    ~H() { mdac_nda_free(p); }
    operator mdac_nda*() const { return p; }
    mdac_nda** out() { return &p; }
};

unsigned env_count(mdac_env* e) {
    unsigned n = 0;
    OK(mdac_env_count(e, &n));
    return n;
}

std::vector<double> coeffs(const mdac_nda* v) {
    size_t n = 0;
    OK(mdac_nda_coeffs(v, nullptr, 0, &n));
    std::vector<double> c(n);
    OK(mdac_nda_coeffs(v, c.data(), n, &n));
    return c;
}

std::vector<double> coeffs(const da::NDA& v) {
    const double* p = v.env_->pool<double>().slot(v.slot_);
    return std::vector<double>(p, p + v.length());
}

// Coefficient equality; trailing zeros do not count.
void require_same(std::vector<double> a, std::vector<double> b) {
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, 0.0);
    b.resize(n, 0.0);
    REQUIRE(a == b);
}

// Test inputs: every coefficient set, constant part c0.
std::vector<double> input(double c0, double scale) {
    std::vector<double> c(da::NDA::full_length());
    c[0] = c0;
    for (size_t i = 1; i < c.size(); ++i) c[i] = scale * ((i % 2) ? 1.0 : -1.0) / double(i + 1);
    return c;
}

mdac_nda* from(mdac_env* env, const std::vector<double>& c) {
    mdac_nda* v = nullptr;
    OK(mdac_nda_from_coeffs(env, c.data(), c.size(), &v));
    return v;
}

} // namespace

// ---- T1.2 ------------------------------------------------------------------

TEST_CASE("capi: mdac_last_error after a failure", "[capi][env]") {
    REQUIRE(mdac_init(2, 2, 0, 0) == MDAC_ERR_VALUE);
    REQUIRE(std::string(mdac_last_error()).size() > 0);
    REQUIRE(std::string(mdac_last_error()) == "poolsize must be positive");
}

TEST_CASE("capi: mdac_init creates the default env and makes it current", "[capi][env]") {
    Envs envs;
    mdac_env* cur = nullptr;
    OK(mdac_env_current(&cur));
    REQUIRE(cur == envs.env);
    unsigned n = 0;
    OK(mdac_env_nvars(envs.env, &n));
    REQUIRE(n == kNvars);
    REQUIRE(env_count(envs.env) == kNvars);  // the base vectors
}

TEST_CASE("capi: mdac_clear retires the default env", "[capi][env]") {
    mdac_env* e = nullptr;
    OK(mdac_init(2, 2, 10, 0));
    OK(mdac_env_default(&e));
    OK(mdac_clear());
    REQUIRE(mdac_env_retired(e) == 1);
    mdac_env* d = nullptr;
    REQUIRE(mdac_env_default(&d) == MDAC_ERR_ENV);
    REQUIRE(std::string(mdac_last_error()).size() > 0);
    OK(mdac_clear());  // no default env: nothing to do
}

TEST_CASE("capi: mdac_env_current fails without a current env", "[capi][env]") {
    Envs envs;
    mdac_env* prev = nullptr;
    OK(mdac_env_exchange(nullptr, &prev));
    mdac_env* cur = nullptr;
    REQUIRE(mdac_env_current(&cur) == MDAC_ERR_ENV);
    REQUIRE(std::string(mdac_last_error()).size() > 0);
    OK(mdac_env_exchange(prev, &prev));
    REQUIRE(prev == nullptr);
    OK(mdac_env_current(&cur));
    REQUIRE(cur == envs.env);
}

TEST_CASE("capi: mdac_env_default", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(2, 2, 10, 0, &e));
    OK(mdac_env_select(e));
    mdac_env* d = nullptr;
    OK(mdac_env_default(&d));
    REQUIRE(d == envs.env);
    OK(mdac_env_close(e));
}

TEST_CASE("capi: mdac_env_make keeps the current env", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(3, 2, 20, 1, &e));
    mdac_env* cur = nullptr;
    OK(mdac_env_current(&cur));
    REQUIRE(cur == envs.env);
    REQUIRE(e != envs.env);
    unsigned n = 0;
    OK(mdac_env_max_order(e, &n));
    REQUIRE(n == 3);
    OK(mdac_env_poolsize(e, &n));
    REQUIRE(n == 20);
    REQUIRE(mdac_env_make(3, 2, 0, 0, &e) == MDAC_ERR_VALUE);
    OK(mdac_env_close(e));
}

TEST_CASE("capi: mdac_env_select", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(2, 2, 10, 0, &e));
    OK(mdac_env_select(e));
    mdac_env* cur = nullptr;
    OK(mdac_env_current(&cur));
    REQUIRE(cur == e);
    OK(mdac_env_select(envs.env));
    OK(mdac_env_close(e));
    REQUIRE(mdac_env_select(e) == MDAC_ERR_ENV);
    REQUIRE(std::string(mdac_last_error()) == "DA environment has been cleared");
}

TEST_CASE("capi: mdac_env_exchange", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(2, 2, 10, 0, &e));
    mdac_env* prev = nullptr;
    OK(mdac_env_exchange(e, &prev));
    REQUIRE(prev == envs.env);
    OK(mdac_env_exchange(prev, &prev));
    REQUIRE(prev == e);
    OK(mdac_env_close(e));
}

TEST_CASE("capi: mdac_env_close", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(2, 2, 10, 0, &e));
    H v;
    OK(mdac_nda_new(e, 1.0, v.out()));
    OK(mdac_env_close(e));
    REQUIRE(mdac_env_retired(e) == 1);
    OK(mdac_env_close(e));  // idempotent
    double x = 0;
    REQUIRE(mdac_nda_con(v, &x) == MDAC_ERR_ENV);
    unsigned n = 0;
    REQUIRE(mdac_env_count(e, &n) == MDAC_ERR_ENV);

    // Closing the default env clears it.
    mdac_env* d = envs.env;
    OK(mdac_env_close(d));
    REQUIRE(mdac_env_retired(d) == 1);
    REQUIRE(mdac_env_default(&d) == MDAC_ERR_ENV);
}

TEST_CASE("capi: env queries", "[capi][env]") {
    Envs envs;
    unsigned n = 0;
    OK(mdac_env_order(envs.env, &n));
    REQUIRE(n == kOrder);
    OK(mdac_env_max_order(envs.env, &n));
    REQUIRE(n == kOrder);
    OK(mdac_env_nvars(envs.env, &n));
    REQUIRE(n == kNvars);
    OK(mdac_env_full_length(envs.env, &n));
    REQUIRE(n == unsigned(da::NDA::full_length()));
    OK(mdac_env_poolsize(envs.env, &n));
    REQUIRE(n == kPool);
    OK(mdac_env_count(envs.env, &n));
    REQUIRE(n == kNvars);
    OK(mdac_env_remain(envs.env, &n));
    REQUIRE(n == kPool - kNvars);
    REQUIRE(mdac_env_retired(envs.env) == 0);
}

TEST_CASE("capi: mdac_env_change_order and mdac_env_restore_order", "[capi][env]") {
    Envs envs;
    int ok = 0;
    OK(mdac_env_change_order(envs.env, 2, &ok));
    REQUIRE(ok == 1);
    unsigned n = 0;
    OK(mdac_env_order(envs.env, &n));
    REQUIRE(n == 2);
    OK(mdac_env_max_order(envs.env, &n));
    REQUIRE(n == kOrder);
    OK(mdac_env_full_length(envs.env, &n));
    REQUIRE(n == 10);  // C(3 + 2, 2)
    OK(mdac_env_change_order(envs.env, kOrder + 1, &ok));
    REQUIRE(ok == 0);
    OK(mdac_env_restore_order(envs.env));
    OK(mdac_env_order(envs.env, &n));
    REQUIRE(n == kOrder);
}

TEST_CASE("capi: mdac_get_eps and mdac_set_eps", "[capi][env]") {
    const double eps = mdac_get_eps();
    REQUIRE(eps > 0);
    OK(mdac_set_eps(1e-10));
    REQUIRE(mdac_get_eps() == 1e-10);
    REQUIRE(mdac_set_eps(0.0) == MDAC_ERR_VALUE);
    REQUIRE(mdac_set_eps(NAN) == MDAC_ERR_VALUE);
    OK(mdac_set_eps(eps));
}

TEST_CASE("capi: pool exhaustion gives MDAC_ERR_POOL", "[capi][env]") {
    Envs envs(kNvars + 2);
    H a, b, c;
    OK(mdac_nda_new(envs.env, 1.0, a.out()));
    OK(mdac_nda_new(envs.env, 2.0, b.out()));
    REQUIRE(mdac_nda_new(envs.env, 3.0, c.out()) == MDAC_ERR_POOL);
    REQUIRE(std::string(mdac_last_error()) == "Pool::assign: Run out of vectors");
    REQUIRE(mdac_nda_add(a, b, c.out()) == MDAC_ERR_POOL);
}

TEST_CASE("capi: vectors of different envs give MDAC_ERR_ENV", "[capi][env]") {
    Envs envs;
    mdac_env* e = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 10, 0, &e));
    H a, b, c;
    OK(mdac_nda_new(envs.env, 1.0, a.out()));
    OK(mdac_nda_new(e, 2.0, b.out()));
    REQUIRE(mdac_nda_add(a, b, c.out()) == MDAC_ERR_ENV);
    REQUIRE(mdac_nda_mul_into(a, a, b) == MDAC_ERR_ENV);
    REQUIRE(mdac_nda_add_d_into(b, a, 1.0) == MDAC_ERR_ENV);
    REQUIRE(mdac_nda_exp_into(b, a) == MDAC_ERR_ENV);
    OK(mdac_env_close(e));
}

// ---- T1.3 ------------------------------------------------------------------

TEST_CASE("capi: mdac_nda_new", "[capi][nda]") {
    Envs envs;
    H v;
    OK(mdac_nda_new(envs.env, 2.5, v.out()));
    require_same(coeffs(v), coeffs(da::NDA(2.5)));
    REQUIRE(env_count(envs.env) == kNvars + 1);
}

TEST_CASE("capi: mdac_nda_from_coeffs", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    REQUIRE(coeffs(v) == c);
    require_same(coeffs(v), coeffs(da::NDA(c)));
    c.push_back(1.0);
    H w;
    REQUIRE(mdac_nda_from_coeffs(envs.env, c.data(), c.size(), w.out()) == MDAC_ERR_VALUE);
}

TEST_CASE("capi: mdac_nda_var", "[capi][nda]") {
    Envs envs;
    for (unsigned i = 0; i < kNvars; ++i) {
        H v;
        OK(mdac_nda_var(envs.env, i, v.out()));
        require_same(coeffs(v), coeffs(da::da_base(i)));
    }
    H v;
    REQUIRE(mdac_nda_var(envs.env, kNvars, v.out()) == MDAC_ERR_INDEX);
}

TEST_CASE("capi: mdac_nda_copy", "[capi][nda]") {
    Envs envs;
    H v(from(envs.env, input(0.3, 1.0))), w;
    OK(mdac_nda_copy(v, w.out()));
    REQUIRE(w.p != v.p);
    REQUIRE(coeffs(w) == coeffs(v));
    REQUIRE(env_count(envs.env) == kNvars + 2);
}

TEST_CASE("capi: mdac_nda_free", "[capi][nda]") {
    Envs envs;
    mdac_nda* v = nullptr;
    OK(mdac_nda_new(envs.env, 1.0, &v));
    REQUIRE(env_count(envs.env) == kNvars + 1);
    mdac_nda_free(v);
    REQUIRE(env_count(envs.env) == kNvars);
    mdac_nda_free(nullptr);
}

TEST_CASE("capi: freeing a vector after mdac_clear is safe, using it fails", "[capi][nda]") {
    mdac_env* e = nullptr;
    OK(mdac_init(kOrder, kNvars, kPool, 0));
    OK(mdac_env_default(&e));
    mdac_nda* v = nullptr;
    OK(mdac_nda_new(e, 1.0, &v));
    OK(mdac_clear());
    double x = 0;
    REQUIRE(mdac_nda_con(v, &x) == MDAC_ERR_ENV);
    mdac_nda* w = nullptr;
    REQUIRE(mdac_nda_exp(v, &w) == MDAC_ERR_ENV);
    REQUIRE(mdac_nda_add_d_into(v, v, 1.0) == MDAC_ERR_ENV);
    REQUIRE(std::string(mdac_last_error()) == "DA environment has been cleared");
    mdac_nda_free(v);
}

TEST_CASE("capi: mdac_nda_env", "[capi][nda]") {
    Envs envs;
    H v;
    OK(mdac_nda_new(envs.env, 1.0, v.out()));
    REQUIRE(mdac_nda_env(v) == envs.env);
}

TEST_CASE("capi: mdac_nda_con and mdac_nda_set_con", "[capi][nda]") {
    Envs envs;
    H v(from(envs.env, input(0.3, 1.0)));
    double x = 0;
    OK(mdac_nda_con(v, &x));
    REQUIRE(x == 0.3);
    OK(mdac_nda_set_con(v, 4.0));
    REQUIRE(coeffs(v) == std::vector<double>{4.0});
}

TEST_CASE("capi: mdac_nda_length, mdac_nda_nterms and mdac_nda_norm", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 3.0);
    c[2] = 0.0;
    H v(from(envs.env, c));
    da::NDA r(c);
    size_t n = 0;
    OK(mdac_nda_length(v, &n));
    REQUIRE(n == r.length());
    OK(mdac_nda_nterms(v, &n));
    REQUIRE(n == size_t(r.n_element()));
    REQUIRE(n == c.size() - 1);
    double x = 0;
    OK(mdac_nda_norm(v, &x));
    REQUIRE(x == r.norm());
}

TEST_CASE("capi: mdac_nda_coeffs size query", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    double buf[4] = {0, 0, 0, -1};
    size_t n = 0;
    OK(mdac_nda_coeffs(v, buf, 3, &n));
    REQUIRE(n == c.size());
    REQUIRE(buf[0] == c[0]);
    REQUIRE(buf[2] == c[2]);
    REQUIRE(buf[3] == -1);  // not written past cap
}

TEST_CASE("capi: mdac_nda_coeff and mdac_nda_set_coeff", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    da::NDA r(c);
    const int exps[] = {1, 0, 2};
    double x = 0;
    OK(mdac_nda_coeff(v, exps, 3, &x));
    REQUIRE(x == r.element(std::vector<int>{1, 0, 2}));
    OK(mdac_nda_set_coeff(v, exps, 3, 7.0));
    OK(mdac_nda_coeff(v, exps, 3, &x));
    REQUIRE(x == 7.0);
    r.set_element(std::vector<int>{1, 0, 2}, 7.0);
    require_same(coeffs(v), coeffs(r));
    const int bad[] = {1, -1, 0};
    REQUIRE(mdac_nda_coeff(v, bad, 3, &x) == MDAC_ERR_VALUE);
    REQUIRE(mdac_nda_set_coeff(v, bad, 3, 1.0) == MDAC_ERR_VALUE);
}

TEST_CASE("capi: mdac_nda_index_term", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    da::NDA r(c);
    for (size_t i = 0; i < c.size(); ++i) {
        int exps[kNvars];
        double x = 0;
        OK(mdac_nda_index_term(v, i, exps, &x));
        std::vector<unsigned> re;
        double rx = 0;
        r.element(unsigned(i), re, rx);
        REQUIRE(x == rx);
        REQUIRE(std::vector<unsigned>(exps, exps + kNvars) == re);
    }
    int exps[kNvars];
    double x = 0;
    REQUIRE(mdac_nda_index_term(v, c.size(), exps, &x) == MDAC_ERR_INDEX);
}

TEST_CASE("capi: mdac_nda_iszero", "[capi][nda]") {
    Envs envs;
    H v;
    OK(mdac_nda_new(envs.env, 1e-8, v.out()));
    int z = -1;
    OK(mdac_nda_iszero(v, 1e-6, &z));
    REQUIRE(z == 1);
    OK(mdac_nda_iszero(v, 1e-10, &z));
    REQUIRE(z == 0);
}

TEST_CASE("capi: mdac_nda_clean", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    da::NDA r(c);
    OK(mdac_nda_clean(v, 0.1));
    r.clean(0.1);
    require_same(coeffs(v), coeffs(r));
}

TEST_CASE("capi: mdac_nda_reset", "[capi][nda]") {
    Envs envs;
    H v(from(envs.env, input(0.3, 1.0)));
    OK(mdac_nda_reset(v));
    for (double x : coeffs(v)) REQUIRE(x == 0.0);
}

TEST_CASE("capi: mdac_nda_to_string", "[capi][nda]") {
    Envs envs;
    std::vector<double> c = input(0.3, 1.0);
    H v(from(envs.env, c));
    std::ostringstream os;
    os << da::NDA(c);
    size_t n = 0;
    OK(mdac_nda_to_string(v, nullptr, 0, &n));
    REQUIRE(n == os.str().size() + 1);
    std::string s(n, 'x');
    OK(mdac_nda_to_string(v, s.data(), n, &n));
    REQUIRE(s.c_str() == os.str());
    char small[5];
    OK(mdac_nda_to_string(v, small, sizeof small, &n));
    REQUIRE(std::string(small) == os.str().substr(0, 4));
}

// ---- T1.4 ------------------------------------------------------------------

namespace {

using Alloc = mdac_status (*)(const mdac_nda*, const mdac_nda*, mdac_nda**);
using AllocD = mdac_status (*)(const mdac_nda*, double, mdac_nda**);
using DAlloc = mdac_status (*)(double, const mdac_nda*, mdac_nda**);
using Into = mdac_status (*)(mdac_nda*, const mdac_nda*, const mdac_nda*);
using IntoD = mdac_status (*)(mdac_nda*, const mdac_nda*, double);
using DInto = mdac_status (*)(mdac_nda*, double, const mdac_nda*);
using Ref = da::NDA (*)(const da::NDA&, const da::NDA&);
using RefD = da::NDA (*)(const da::NDA&, double);
using DRef = da::NDA (*)(double, const da::NDA&);

struct BinOp {
    const char* name;
    Alloc f;
    AllocD fd;
    DAlloc df;
    Into fi;
    IntoD fdi;
    DInto dfi;
    Ref r;
    RefD rd;
    DRef dr;
};

#define BINOP(name, op)                                                              \
    BinOp{#name, mdac_nda_##name, mdac_nda_##name##_d, mdac_nda_d##name,             \
          mdac_nda_##name##_into, mdac_nda_##name##_d_into, mdac_nda_d##name##_into, \
          [](const da::NDA& a, const da::NDA& b) { return a op b; },                 \
          [](const da::NDA& a, double x) { return a op x; },                         \
          [](double x, const da::NDA& a) { return x op a; }}

const BinOp kBinOps[] = {BINOP(add, +), BINOP(sub, -), BINOP(mul, *), BINOP(div, /)};

} // namespace

TEST_CASE("capi: NDA arithmetic equals C++", "[capi][arith]") {
    Envs envs;
    const std::vector<double> ca = input(0.3, 1.0), cb = input(1.7, -0.5);
    H a(from(envs.env, ca)), b(from(envs.env, cb));
    const da::NDA ra(const_cast<std::vector<double>&>(ca)), rb(const_cast<std::vector<double>&>(cb));
    const double x = 1.25;

    for (const BinOp& op : kBinOps) {
        INFO(op.name);
        {
            H o;
            OK(op.f(a, b, o.out()));
            require_same(coeffs(o), coeffs(op.r(ra, rb)));
            H od;
            OK(op.fd(a, x, od.out()));
            require_same(coeffs(od), coeffs(op.rd(ra, x)));
            H dO;
            OK(op.df(x, a, dO.out()));
            require_same(coeffs(dO), coeffs(op.dr(x, ra)));
        }

        // In place into a vector holding other values; no new slot.
        H o(from(envs.env, input(-2.0, 3.0)));
        const unsigned n = env_count(envs.env);
        OK(op.fi(o, a, b));
        REQUIRE(env_count(envs.env) == n);
        require_same(coeffs(o), coeffs(op.r(ra, rb)));
        OK(op.fdi(o, a, x));
        REQUIRE(env_count(envs.env) == n);
        require_same(coeffs(o), coeffs(op.rd(ra, x)));
        OK(op.dfi(o, x, a));
        REQUIRE(env_count(envs.env) == n);
        require_same(coeffs(o), coeffs(op.dr(x, ra)));

        // Aliasing: out == a, out == b, out == a == b.
        H a1(from(envs.env, ca)), b1(from(envs.env, cb));
        OK(op.fi(a1, a1, b));
        require_same(coeffs(a1), coeffs(op.r(ra, rb)));
        OK(op.fi(b1, a, b1));
        require_same(coeffs(b1), coeffs(op.r(ra, rb)));
        H s(from(envs.env, ca));
        OK(op.fi(s, s, s));
        require_same(coeffs(s), coeffs(op.r(ra, ra)));
        H s1(from(envs.env, ca)), s2(from(envs.env, ca));
        OK(op.fdi(s1, s1, x));
        require_same(coeffs(s1), coeffs(op.rd(ra, x)));
        OK(op.dfi(s2, x, s2));
        require_same(coeffs(s2), coeffs(op.dr(x, ra)));
        REQUIRE(env_count(envs.env) == n + 5);
    }
}

TEST_CASE("capi: division by zero gives MDAC_ERR_VALUE", "[capi][arith]") {
    Envs envs;
    H a(from(envs.env, input(0.3, 1.0))), o;
    REQUIRE(mdac_nda_div_d(a, 0.0, o.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_nda_div_d_into(a, a, 0.0) == MDAC_ERR_VALUE);
}

TEST_CASE("capi: neg, pow_i and pow_d equal C++", "[capi][arith]") {
    Envs envs;
    const std::vector<double> ca = input(1.3, 1.0);
    H a(from(envs.env, ca));
    const da::NDA ra(const_cast<std::vector<double>&>(ca));
    H o(from(envs.env, input(-2.0, 3.0)));
    const unsigned n = env_count(envs.env);

    H neg, pi, pd;
    OK(mdac_nda_neg(a, neg.out()));
    require_same(coeffs(neg), coeffs(-ra));
    OK(mdac_nda_neg_into(o, a));
    require_same(coeffs(o), coeffs(-ra));
    OK(mdac_nda_pow_i(a, 3, pi.out()));
    require_same(coeffs(pi), coeffs(da::pow(ra, 3)));
    OK(mdac_nda_pow_i_into(o, a, 3));
    require_same(coeffs(o), coeffs(da::pow(ra, 3)));
    OK(mdac_nda_pow_d(a, 0.3, pd.out()));
    require_same(coeffs(pd), coeffs(da::pow(ra, 0.3)));
    OK(mdac_nda_pow_d_into(o, a, 0.3));
    require_same(coeffs(o), coeffs(da::pow(ra, 0.3)));
    REQUIRE(env_count(envs.env) == n + 3);

    H s(from(envs.env, ca));
    OK(mdac_nda_pow_i_into(s, s, 2));
    require_same(coeffs(s), coeffs(da::pow(ra, 2)));
    OK(mdac_nda_neg_into(s, s));
    require_same(coeffs(s), coeffs(-da::pow(ra, 2)));
}

// ---- T1.5 ------------------------------------------------------------------

namespace {

struct Func {
    const char* name;
    mdac_status (*f)(const mdac_nda*, mdac_nda**);
    mdac_status (*fi)(mdac_nda*, const mdac_nda*);
    da::NDA (*r)(const da::NDA&);
    double c0;  // a constant part inside the domain
};

#define FUNC(name, c0) Func{#name, mdac_nda_##name, mdac_nda_##name##_into, da::name, c0}

const Func kFuncs[] = {
    FUNC(sqrt, 1.3),  FUNC(exp, 0.3),   FUNC(log, 1.3),   FUNC(sin, 0.3),
    FUNC(cos, 0.3),   FUNC(tan, 0.3),   FUNC(asin, 0.3),  FUNC(acos, 0.3),
    FUNC(atan, 0.3),  FUNC(sinh, 0.3),  FUNC(cosh, 0.3),  FUNC(tanh, 0.3),
    FUNC(asinh, 0.3), FUNC(acosh, 1.7), FUNC(atanh, 0.3), FUNC(erf, 0.3),
};

} // namespace

TEST_CASE("capi: NDA math functions equal C++", "[capi][func]") {
    Envs envs;
    for (const Func& fn : kFuncs) {
        INFO(fn.name);
        const std::vector<double> ca = input(fn.c0, 0.2);
        H a(from(envs.env, ca)), o;
        const da::NDA ra(const_cast<std::vector<double>&>(ca));
        const std::vector<double> want = coeffs(fn.r(ra));
        OK(fn.f(a, o.out()));
        require_same(coeffs(o), want);

        H into(from(envs.env, input(-2.0, 3.0)));
        const unsigned n = env_count(envs.env);
        OK(fn.fi(into, a));
        REQUIRE(env_count(envs.env) == n);
        require_same(coeffs(into), want);
        OK(fn.fi(a, a));
        require_same(coeffs(a), want);
    }
}

TEST_CASE("capi: mdac_nda_abs", "[capi][func]") {
    Envs envs;
    const std::vector<double> ca = input(0.3, -4.0);
    H a(from(envs.env, ca));
    double x = 0;
    OK(mdac_nda_abs(a, &x));
    REQUIRE(x == da::abs(da::NDA(const_cast<std::vector<double>&>(ca))));
}

TEST_CASE("capi: asin outside its domain gives MDAC_ERR_VALUE", "[capi][func]") {
    Envs envs;
    H a, o;
    OK(mdac_nda_new(envs.env, 2.0, a.out()));
    REQUIRE(mdac_nda_asin(a, o.out()) == MDAC_ERR_VALUE);
    REQUIRE(o.p == nullptr);
    REQUIRE(mdac_nda_asin_into(a, a) == MDAC_ERR_VALUE);
    REQUIRE(std::string(mdac_last_error()).find("asin") != std::string::npos);
}

// ---- T3.1 ------------------------------------------------------------------

namespace {

// Owns a C API list.
struct L {
    mdac_ndalist* p = nullptr;
    L() = default;
    L(const L&) = delete;
    L& operator=(const L&) = delete;
    ~L() { mdac_ndalist_free(p); }
    operator mdac_ndalist*() const { return p; }
    mdac_ndalist** out() { return &p; }
};

std::vector<double> elem(const mdac_ndalist* l, size_t i) {
    H v;
    OK(mdac_ndalist_get(l, i, v.out()));
    return coeffs(v);
}

// x = 1 + x0 + 2 x1 + 5 x2 as in python/tests/test_algorithms.py, and its C++ twin.
struct X {
    H h;
    da::NDA r;
    explicit X(mdac_env* e) : r(1.0 + da::da_base(0) + 2.0 * da::da_base(1) + 5.0 * da::da_base(2)) {
        h.p = from(e, coeffs(r));
    }
};

// A C API list of copies of hs.
void make_list(std::initializer_list<const mdac_nda*> hs, L& l) {
    std::vector<const mdac_nda*> v(hs);
    OK(mdac_ndalist_from(v.data(), v.size(), l.out()));
}

} // namespace

TEST_CASE("capi: NDA lists", "[capi][list]") {
    Envs envs;
    X x(envs.env);
    H e;
    OK(mdac_nda_exp(x.h, e.out()));
    L l;
    OK(mdac_ndalist_new(l.out()));
    REQUIRE(mdac_ndalist_length(l) == 0);
    REQUIRE(mdac_ndalist_env(l) == nullptr);
    OK(mdac_ndalist_push(l, x.h));
    OK(mdac_ndalist_push(l, e));
    REQUIRE(mdac_ndalist_length(l) == 2);
    REQUIRE(mdac_ndalist_env(l) == envs.env);
    require_same(elem(l, 1), coeffs(da::exp(x.r)));

    const unsigned n = env_count(envs.env);
    OK(mdac_ndalist_set(l, 0, e));                  // copies into the element's slot
    REQUIRE(env_count(envs.env) == n);
    require_same(elem(l, 0), coeffs(da::exp(x.r)));

    H v;
    REQUIRE(mdac_ndalist_get(l, 2, v.out()) == MDAC_ERR_INDEX);
    REQUIRE(mdac_ndalist_set(l, 2, e) == MDAC_ERR_INDEX);

    L m;
    make_list({x.h, e, x.h}, m);
    REQUIRE(mdac_ndalist_length(m) == 3);
    require_same(elem(m, 2), coeffs(x.r));
    REQUIRE(env_count(envs.env) == n + 3);
    mdac_ndalist_free(m.p);
    m.p = nullptr;
    REQUIRE(env_count(envs.env) == n);
    mdac_ndalist_free(nullptr);
}

TEST_CASE("capi: a list holds vectors of one env", "[capi][list]") {
    Envs envs;
    mdac_env* other = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 50, 0, &other));
    H a, b;
    OK(mdac_nda_new(envs.env, 1.0, a.out()));
    OK(mdac_nda_new(other, 2.0, b.out()));
    L l;
    make_list({a}, l);
    REQUIRE(mdac_ndalist_push(l, b) == MDAC_ERR_ENV);
    REQUIRE(mdac_ndalist_set(l, 0, b) == MDAC_ERR_ENV);
    std::vector<const mdac_nda*> ab{a, b};
    L m;
    REQUIRE(mdac_ndalist_from(ab.data(), 2, m.out()) == MDAC_ERR_ENV);
    mdac_nda_free(b.p);
    b.p = nullptr;
    OK(mdac_env_close(other));
}

TEST_CASE("capi: a list outliving its env is freed safely", "[capi][list]") {
    L l;
    {
        OK(mdac_init(kOrder, kNvars, kPool, 0));
        mdac_env* e = nullptr;
        OK(mdac_env_default(&e));
        H a;
        OK(mdac_nda_new(e, 1.0, a.out()));
        make_list({a, a}, l);
        OK(mdac_clear());
    }
    H v;
    REQUIRE(mdac_ndalist_get(l, 0, v.out()) == MDAC_ERR_ENV);
}

TEST_CASE("capi: der and integ equal C++", "[capi][algo]") {
    Envs envs;
    X x(envs.env);
    H y, d, i;
    OK(mdac_nda_log(x.h, y.out()));
    const da::NDA ry = da::log(x.r);
    OK(mdac_nda_der(y, 1, d.out()));
    require_same(coeffs(d), coeffs(da::da_der(ry, 1)));
    OK(mdac_nda_integ(y, 1, i.out()));
    require_same(coeffs(i), coeffs(da::da_int(ry, 1)));
    H bad;
    REQUIRE(mdac_nda_der(y, 3, bad.out()) == MDAC_ERR_INDEX);
    REQUIRE(mdac_nda_integ(y, 3, bad.out()) == MDAC_ERR_INDEX);
}

TEST_CASE("capi: substitute equals C++", "[capi][algo]") {
    Envs envs;
    X x(envs.env);
    H e, s, c;
    OK(mdac_nda_exp(x.h, e.out()));
    OK(mdac_nda_sin(x.h, s.out()));
    OK(mdac_nda_cos(x.h, c.out()));
    const da::NDA re = da::exp(x.r);

    H o1, o2, o3;
    OK(mdac_nda_substitute_d(e, 0, 1.0, o1.out()));
    da::NDA r1;
    da::da_substitute_const(re, 0, 1.0, r1);
    require_same(coeffs(o1), coeffs(r1));

    OK(mdac_nda_substitute(e, 0, x.h, o2.out()));
    da::NDA r2;
    da::da_substitute(re, 0, x.r, r2);
    require_same(coeffs(o2), coeffs(r2));

    L sc;
    make_list({s, c}, sc);
    const unsigned ids[] = {0, 1};
    OK(mdac_nda_substitute_multi(e, ids, 2, sc, o3.out()));
    std::vector<unsigned> rids{0, 1};
    std::vector<da::NDA> rsc{da::sin(x.r), da::cos(x.r)};
    da::NDA r3;
    da::da_substitute(re, rids, rsc, r3);
    require_same(coeffs(o3), coeffs(r3));

    H bad;
    const unsigned dup[] = {0, 0}, far[] = {0, 3};
    REQUIRE(mdac_nda_substitute_multi(e, dup, 2, sc, bad.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_nda_substitute_multi(e, far, 2, sc, bad.out()) == MDAC_ERR_INDEX);
    REQUIRE(mdac_nda_substitute_multi(e, ids, 1, sc, bad.out()) == MDAC_ERR_VALUE);

    L m, out;
    H sh;
    OK(mdac_nda_sinh(x.h, sh.out()));
    make_list({x.h, e, sh}, m);
    OK(mdac_ndalist_substitute(m, ids, 2, sc, out.out()));
    std::vector<da::NDA> rm{x.r, re, da::sinh(x.r)}, rout(3);
    da::da_substitute(rm, rids, rsc, rout);
    REQUIRE(mdac_ndalist_length(out) == 3);
    for (size_t k = 0; k < 3; ++k) require_same(elem(out, k), coeffs(rout[k]));
}

TEST_CASE("capi: compose equals C++", "[capi][algo]") {
    Envs envs;
    X x(envs.env);
    H e, sh, s, c, t;
    OK(mdac_nda_exp(x.h, e.out()));
    OK(mdac_nda_sinh(x.h, sh.out()));
    OK(mdac_nda_sin(x.h, s.out()));
    OK(mdac_nda_cos(x.h, c.out()));
    OK(mdac_nda_tan(x.h, t.out()));
    L m, v, out;
    make_list({x.h, e, sh}, m);
    make_list({s, c, t}, v);
    std::vector<da::NDA> rm{x.r, da::exp(x.r), da::sinh(x.r)};
    std::vector<da::NDA> rv{da::sin(x.r), da::cos(x.r), da::tan(x.r)}, rout(3);

    OK(mdac_ndalist_compose(m, v, out.out()));
    da::da_composition(rm, rv, rout);
    for (size_t k = 0; k < 3; ++k) require_same(elem(out, k), coeffs(rout[k]));

    std::vector<double> pt{0.1, -0.2, 0.3}, got(3), want(3);
    OK(mdac_ndalist_compose_d(m, pt.data(), 3, got.data()));
    da::da_composition(rm, pt, want);
    REQUIRE(got == want);

    std::vector<std::complex<double>> zpt{{0.1, 0.2}, {0.0, -0.1}, {0.3, 0.0}}, zgot(3), zwant(3);
    OK(mdac_ndalist_compose_z(m, reinterpret_cast<const double*>(zpt.data()), 3,
                              reinterpret_cast<double*>(zgot.data())));
    da::da_composition(rm, zpt, zwant);
    REQUIRE(zgot == zwant);

    L bad;
    REQUIRE(mdac_ndalist_compose_d(m, pt.data(), 2, got.data()) == MDAC_ERR_VALUE);
    L two;
    make_list({s, c}, two);
    REQUIRE(mdac_ndalist_compose(m, two, bad.out()) == MDAC_ERR_VALUE);
}

TEST_CASE("capi: inv_map equals C++ and rejects bad maps", "[capi][algo]") {
    OK(mdac_init(5, 2, 1000, 0));
    da::da_init(5, 2, 1000);
    mdac_env* env = nullptr;
    OK(mdac_env_default(&env));
    {
        const da::NDA x = da::da_base(0), y = da::da_base(1);
        std::vector<da::NDA> rm{2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y,
                                -0.4 * x + 1.5 * y + 0.2 * y * y}, rinv(2);
        da::inv_map(rm, 2, rinv);
        H a(from(env, coeffs(rm[0]))), b(from(env, coeffs(rm[1])));
        L m, inv;
        make_list({a, b}, m);
        OK(mdac_ndalist_inv_map(m, 2, inv.out()));
        REQUIRE(mdac_ndalist_length(inv) == 2);
        for (size_t k = 0; k < 2; ++k) require_same(elem(inv, k), coeffs(rinv[k]));

        L bad;
        REQUIRE(mdac_ndalist_inv_map(m, 3, bad.out()) == MDAC_ERR_VALUE);
        H a1;
        OK(mdac_nda_add_d(a, 1.0, a1.out()));
        L shifted;
        make_list({a1, b}, shifted);
        REQUIRE(mdac_ndalist_inv_map(shifted, 2, bad.out()) == MDAC_ERR_VALUE);
    }
    mdac_clear();
    da::da_clear();
}

TEST_CASE("capi: evaluate_map and eval equal compose", "[capi][algo]") {
    Envs envs;
    X x(envs.env);
    H e, sh;
    OK(mdac_nda_exp(x.h, e.out()));
    OK(mdac_nda_sinh(x.h, sh.out()));
    L m;
    make_list({x.h, e, sh}, m);
    const size_t np = 300;
    std::vector<double> pts(np * 3), out(np * 3);
    for (size_t i = 0; i < pts.size(); ++i) pts[i] = std::sin(0.37 * double(i));
    OK(mdac_ndalist_evaluate_map(m, pts.data(), np, out.data()));
    for (size_t p = 0; p < np; ++p) {
        double want[3];
        OK(mdac_ndalist_compose_d(m, pts.data() + 3 * p, 3, want));
        for (size_t c = 0; c < 3; ++c) REQUIRE(out[3 * p + c] == want[c]);
        double v = 0;
        OK(mdac_nda_eval(e, pts.data() + 3 * p, 3, &v));
        REQUIRE(v == want[1]);
    }
    double v = 0;
    REQUIRE(mdac_nda_eval(e, pts.data(), 2, &v) == MDAC_ERR_VALUE);
}

TEST_CASE("capi: mdac_exponents", "[capi][algo]") {
    Envs envs;
    size_t n = 0;
    OK(mdac_exponents(envs.env, nullptr, 0, &n));
    const size_t len = da::NDA::full_length();
    REQUIRE(n == len * kNvars);
    std::vector<int> buf(n);
    OK(mdac_exponents(envs.env, buf.data(), n, &n));
    int exps[kNvars];
    for (size_t i = 0; i < len; ++i) {
        H v(from(envs.env, std::vector<double>(len, 1.0)));
        double c = 0;
        OK(mdac_nda_index_term(v, i, exps, &c));
        for (size_t j = 0; j < kNvars; ++j) REQUIRE(buf[i * kNvars + j] == exps[j]);
    }
}

// ---- T4.1 ------------------------------------------------------------------

namespace {

using da::get_imag;
using da::get_real;

// Owns a C API CNDA.
struct C {
    mdac_cnda* p = nullptr;
    C() = default;
    explicit C(mdac_cnda* q) : p(q) {}
    C(const C&) = delete;
    C& operator=(const C&) = delete;
    ~C() { mdac_cnda_free(p); }
    operator mdac_cnda*() const { return p; }
    mdac_cnda** out() { return &p; }
};

// Owns a C API CNDA list.
struct CL {
    mdac_cndalist* p = nullptr;
    CL() = default;
    CL(const CL&) = delete;
    CL& operator=(const CL&) = delete;
    ~CL() { mdac_cndalist_free(p); }
    operator mdac_cndalist*() const { return p; }
    mdac_cndalist** out() { return &p; }
};

// A C API copy of the C++ CNDA r.
mdac_cnda* cfrom(mdac_env* e, const da::CNDA& r) {
    H re(from(e, coeffs(get_real(r)))), im(from(e, coeffs(get_imag(r))));
    mdac_cnda* c = nullptr;
    OK(mdac_cnda_new(re, im, &c));
    return c;
}

void require_same(const mdac_cnda* c, const da::CNDA& r) {
    H re, im;
    OK(mdac_cnda_real(c, re.out()));
    OK(mdac_cnda_imag(c, im.out()));
    require_same(coeffs(re), coeffs(get_real(r)));
    require_same(coeffs(im), coeffs(get_imag(r)));
}

// The allocating form and the _into form (into o, whose slots it keeps) both give want.
template <class F, class G>
void both(mdac_env* e, mdac_cnda* o, F alloc, G into, const da::CNDA& want) {
    {
        C r;
        OK(alloc(r.out()));
        require_same(r, want);
    }
    const unsigned n = env_count(e);
    OK(into(o));
    REQUIRE(env_count(e) == n);
    require_same(o, want);
}

da::NDA nda(double c0, double scale) {
    std::vector<double> c = input(c0, scale);
    return da::NDA(c);
}

#define ALLOC(expr) [&](mdac_cnda** r_) { return expr; }
#define INTO(expr) [&](mdac_cnda* t_) { return expr; }

#define CNDA_OP(name, op)                                                                        \
    SECTION(#name) {                                                                             \
        both(e, o, ALLOC(mdac_cnda_##name(a, b, r_)), INTO(mdac_cnda_##name##_into(t_, a, b)), ra op rb); \
        both(e, o, ALLOC(mdac_cnda_##name##_n(a, x, r_)), INTO(mdac_cnda_##name##_n_into(t_, a, x)), ra op rx); \
        both(e, o, ALLOC(mdac_cnda_n##name(x, a, r_)), INTO(mdac_cnda_n##name##_into(t_, x, a)), rx op ra); \
        both(e, o, ALLOC(mdac_cnda_##name##_d(a, d, r_)), INTO(mdac_cnda_##name##_d_into(t_, a, d)), ra op d); \
        both(e, o, ALLOC(mdac_cnda_d##name(d, a, r_)), INTO(mdac_cnda_d##name##_into(t_, d, a)), d op ra); \
        both(e, o, ALLOC(mdac_cnda_##name##_z(a, z.real(), z.imag(), r_)),                        \
             INTO(mdac_cnda_##name##_z_into(t_, a, z.real(), z.imag())), ra op z);               \
        both(e, o, ALLOC(mdac_cnda_z##name(z.real(), z.imag(), a, r_)),                           \
             INTO(mdac_cnda_z##name##_into(t_, z.real(), z.imag(), a)), z op ra);                \
        both(e, o, ALLOC(mdac_nda_##name##_z(x, z.real(), z.imag(), r_)),                         \
             INTO(mdac_nda_##name##_z_into(t_, x, z.real(), z.imag())), rx op z);                \
        both(e, o, ALLOC(mdac_nda_z##name(z.real(), z.imag(), x, r_)),                            \
             INTO(mdac_nda_z##name##_into(t_, z.real(), z.imag(), x)), z op rx);                 \
        C a1(cfrom(e, ra)), b1(cfrom(e, rb));                                                    \
        OK(mdac_cnda_##name##_into(a1, a1, b));                                                  \
        require_same(a1, ra op rb);                                                              \
        OK(mdac_cnda_##name##_into(b1, a, b1));                                                  \
        require_same(b1, ra op rb);                                                              \
    }

} // namespace

TEST_CASE("capi: CNDA lifecycle and parts", "[capi][cnda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::NDA rre = nda(0.3, 1.0), rim = nda(1.7, -0.5);
    H re(from(e, coeffs(rre))), im(from(e, coeffs(rim)));

    C c, r0, z;
    OK(mdac_cnda_new(re, im, c.out()));
    require_same(c, da::CNDA(rre, rim));
    OK(mdac_cnda_new(re, nullptr, r0.out()));
    require_same(r0, da::CNDA(rre, da::NDA(0.0)));
    OK(mdac_cnda_new_z(e, 1.5, -2.0, z.out()));
    require_same(z, da::CNDA(da::NDA(1.5), da::NDA(-2.0)));
    REQUIRE(mdac_cnda_env(c) == e);

    C cp;
    OK(mdac_cnda_copy(c, cp.out()));
    require_same(cp, da::CNDA(rre, rim));

    // The parts are copies; the setters copy into the CNDA's own slots.
    H p;
    OK(mdac_cnda_real(c, p.out()));
    OK(mdac_nda_set_con(p, 9.0));
    require_same(c, da::CNDA(rre, rim));
    const unsigned n = env_count(e);
    OK(mdac_cnda_set_real(c, im));
    OK(mdac_cnda_set_imag(c, re));
    REQUIRE(env_count(e) == n);
    require_same(c, da::CNDA(rim, rre));

    size_t len = 0;
    OK(mdac_cnda_to_string(c, nullptr, 0, &len));
    std::string s(len, '\0');
    OK(mdac_cnda_to_string(c, s.data(), len, &len));
    std::ostringstream os;
    os << da::CNDA(rim, rre);
    // The first line names the slots, which differ between the two envs.
    const std::string want = os.str(), got = s.c_str();
    REQUIRE(got.find("V [") != std::string::npos);
    REQUIRE(got.substr(got.find('\n')) == want.substr(want.find('\n')));
}

TEST_CASE("capi: CNDA parts and operands must share the env", "[capi][cnda]") {
    Envs envs;
    mdac_env* other = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 50, 0, &other));
    H a(from(envs.env, input(0.3, 1.0))), b;
    OK(mdac_nda_new(other, 1.0, b.out()));
    C c, o;
    REQUIRE(mdac_cnda_new(a, b, c.out()) == MDAC_ERR_ENV);
    OK(mdac_cnda_new(a, nullptr, c.out()));
    REQUIRE(mdac_cnda_set_imag(c, b) == MDAC_ERR_ENV);
    REQUIRE(mdac_cnda_add_n(c, b, o.out()) == MDAC_ERR_ENV);
    REQUIRE(mdac_cnda_mul_n_into(c, c, b) == MDAC_ERR_ENV);
    OK(mdac_env_close(other));
}

TEST_CASE("capi: freeing a CNDA after mdac_clear is safe, using it fails", "[capi][cnda]") {
    mdac_cnda* c = nullptr;
    {
        Envs envs;
        OK(mdac_cnda_new_z(envs.env, 1.0, 2.0, &c));
    }
    C o;
    REQUIRE(mdac_cnda_exp(c, o.out()) == MDAC_ERR_ENV);
    mdac_cnda_free(c);
}

TEST_CASE("capi: CNDA arithmetic equals C++", "[capi][cnda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::CNDA ra(nda(0.3, 1.0), nda(1.7, -0.5)), rb(nda(1.1, 0.5), nda(-0.4, 0.25));
    const da::NDA rx = nda(0.9, -1.0);
    const double d = 1.25;
    const std::complex<double> z(0.5, -2.0);
    C a(cfrom(e, ra)), b(cfrom(e, rb)), o(cfrom(e, da::CNDA(nda(-2.0, 3.0), nda(1.0, 1.0))));
    H x(from(e, coeffs(rx)));

    CNDA_OP(add, +)
    CNDA_OP(sub, -)
    CNDA_OP(mul, *)
    CNDA_OP(div, /)
    SECTION("neg and pow") {
        both(e, o, ALLOC(mdac_cnda_neg(a, r_)), INTO(mdac_cnda_neg_into(t_, a)), -ra);
        both(e, o, ALLOC(mdac_cnda_pow_i(a, 3, r_)), INTO(mdac_cnda_pow_i_into(t_, a, 3)), da::pow(ra, 3));
        both(e, o, ALLOC(mdac_cnda_pow_d(a, 0.3, r_)), INTO(mdac_cnda_pow_d_into(t_, a, 0.3)),
             da::pow(ra, 0.3));
        OK(mdac_cnda_neg_into(a, a));
        require_same(a, -ra);
    }
    SECTION("division by zero") {
        C r;
        REQUIRE(mdac_cnda_div_d(a, 0.0, r.out()) == MDAC_ERR_VALUE);
        REQUIRE(mdac_cnda_div_d_into(o, a, 0.0) == MDAC_ERR_VALUE);
    }
}

TEST_CASE("capi: CNDA math functions equal C++", "[capi][cnda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::CNDA ra(nda(0.3, 0.2), nda(0.2, 0.1));
    C a(cfrom(e, ra)), o(cfrom(e, da::CNDA(nda(-2.0, 3.0), nda(1.0, 1.0))));

#define CNDA_FUNC(f) \
    both(e, o, ALLOC(mdac_cnda_##f(a, r_)), INTO(mdac_cnda_##f##_into(t_, a)), da::f(ra));
    CNDA_FUNC(sqrt)
    CNDA_FUNC(exp)
    CNDA_FUNC(log)
    CNDA_FUNC(asin)
    CNDA_FUNC(acos)
    CNDA_FUNC(atan)
    CNDA_FUNC(asinh)
    CNDA_FUNC(acosh)
    CNDA_FUNC(atanh)
#undef CNDA_FUNC

    double x = 0;
    OK(mdac_cnda_abs(a, &x));
    REQUIRE(x == da::abs(ra));
    OK(mdac_cnda_exp_into(a, a));
    require_same(a, da::exp(ra));
}

TEST_CASE("capi: CNDA lists", "[capi][cnda][list]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::CNDA r1(nda(0.3, 1.0), nda(1.7, -0.5)), r2(nda(1.1, 0.5), nda(-0.4, 0.25));
    C c1(cfrom(e, r1)), c2(cfrom(e, r2));

    CL l, empty;
    OK(mdac_cndalist_new(empty.out()));
    REQUIRE(mdac_cndalist_length(empty) == 0);
    REQUIRE(mdac_cndalist_env(empty) == nullptr);
    std::vector<const mdac_cnda*> vs{c1};
    OK(mdac_cndalist_from(vs.data(), vs.size(), l.out()));
    OK(mdac_cndalist_push(l, c2));
    REQUIRE(mdac_cndalist_length(l) == 2);
    REQUIRE(mdac_cndalist_env(l) == e);
    {
        C g;
        OK(mdac_cndalist_get(l, 1, g.out()));
        require_same(g, r2);
    }
    const unsigned n = env_count(e);
    OK(mdac_cndalist_set(l, 1, c1));
    REQUIRE(env_count(e) == n);
    C g;
    OK(mdac_cndalist_get(l, 1, g.out()));
    require_same(g, r1);
    REQUIRE(mdac_cndalist_get(l, 2, g.out()) == MDAC_ERR_INDEX);
    REQUIRE(mdac_cndalist_set(l, 2, c1) == MDAC_ERR_INDEX);

    mdac_env* other = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 50, 0, &other));
    {
        C c3;
        OK(mdac_cnda_new_z(other, 1.0, 0.0, c3.out()));
        REQUIRE(mdac_cndalist_push(l, c3) == MDAC_ERR_ENV);
        REQUIRE(mdac_cndalist_set(l, 0, c3) == MDAC_ERR_ENV);
        vs.push_back(c3);
        CL bad;
        REQUIRE(mdac_cndalist_from(vs.data(), vs.size(), bad.out()) == MDAC_ERR_ENV);
    }
    OK(mdac_env_close(other));
}

TEST_CASE("capi: cd_composition equals C++", "[capi][cnda][algo]") {
    Envs envs;
    mdac_env* e = envs.env;
    // The vectors of test_numeric.cc's "CD FUNCTIONS (numeric)".
    da::NDA t = da::da_base(0) + 2.0 * da::da_base(1) + 3.0 * da::da_base(2);
    da::NDA s = 0.5 * da::da_base(0) + 4.0 * da::da_base(1) + 2.7 * da::da_base(2);
    const da::NDA x1 = da::cos(t), x2 = da::sin(t), x3 = da::cos(s), x4 = da::sin(s);
    const da::CNDA y1(x1, x2), y2(x3, x4);
    std::vector<da::NDA> mmap{x1, x2};
    std::vector<da::CNDA> cnmap{y1, y2, y1 * y2};
    std::vector<da::CNDA> cmmap{da::CNDA(x1, da::exp(x1)), da::CNDA(x2, da::exp(x2))};
    std::vector<da::NDA> nmap{x1, x2, x1 + 0.33 * x2};

    H h1(from(e, coeffs(x1))), h2(from(e, coeffs(x2))), h3(from(e, coeffs(nmap[2])));
    C cy1(cfrom(e, cnmap[0])), cy2(cfrom(e, cnmap[1])), cy3(cfrom(e, cnmap[2]));
    C cm1(cfrom(e, cmmap[0])), cm2(cfrom(e, cmmap[1]));
    L m, n;
    make_list({h1, h2}, m);
    make_list({h1, h2, h3}, n);
    CL cn, cm;
    std::vector<const mdac_cnda*> vcn{cy1, cy2, cy3}, vcm{cm1, cm2};
    OK(mdac_cndalist_from(vcn.data(), vcn.size(), cn.out()));
    OK(mdac_cndalist_from(vcm.data(), vcm.size(), cm.out()));

    auto check = [&](mdac_cndalist* got, const std::vector<da::CNDA>& want) {
        REQUIRE(mdac_cndalist_length(got) == want.size());
        for (size_t k = 0; k < want.size(); ++k) {
            C g;
            OK(mdac_cndalist_get(got, k, g.out()));
            require_same(g, want[k]);
        }
    };
    std::vector<da::CNDA> want(2);
    {
        CL out;
        OK(mdac_ndalist_compose_c(m, cn, out.out()));
        da::cd_composition(mmap, cnmap, want);
        check(out, want);
    }
    {
        CL out;
        OK(mdac_cndalist_compose(cm, cn, out.out()));
        da::cd_composition(cmmap, cnmap, want);
        check(out, want);
    }
    {
        CL out;
        OK(mdac_cndalist_compose_n(cm, n, out.out()));
        da::cd_composition(cmmap, nmap, want);
        check(out, want);
    }

    // An NDA map at complex points equals cd_composition at constant CNDA arguments.
    std::vector<std::complex<double>> pt{{0.1, 0.2}, {-0.3, 0.05}, {0.2, -0.1}}, vals(2);
    OK(mdac_ndalist_compose_z(m, reinterpret_cast<const double*>(pt.data()), 3,
                              reinterpret_cast<double*>(vals.data())));
    std::vector<da::CNDA> cpt;
    for (auto p : pt) cpt.emplace_back(da::NDA(p.real()), da::NDA(p.imag()));
    da::cd_composition(mmap, cpt, want);
    for (size_t k = 0; k < 2; ++k) {
        REQUIRE(vals[k].real() == Approx(get_real(want[k]).con()).epsilon(1e-14).margin(1e-15));
        REQUIRE(vals[k].imag() == Approx(get_imag(want[k]).con()).epsilon(1e-14).margin(1e-15));
    }

    CL bad, two;
    std::vector<const mdac_cnda*> v2{cy1, cy2};
    OK(mdac_cndalist_from(v2.data(), v2.size(), two.out()));
    REQUIRE(mdac_ndalist_compose_c(m, two, bad.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_cndalist_compose_n(cm, m, bad.out()) == MDAC_ERR_VALUE);
}

// ---- Multi-env (T7.1) -----------------------------------------------------------

TEST_CASE("capi: import copies NDA and CNDA into another env", "[capi][multienv]") {
    Envs envs;
    mdac_env *b = nullptr, *small = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 20, 0, &b));
    OK(mdac_env_make(2, 2, 20, 0, &small));
    H x(from(envs.env, input(1.5, 0.7)));
    H y;
    OK(mdac_nda_import(b, x, y.out()));
    REQUIRE(mdac_nda_env(y) == b);
    REQUIRE(env_count(b) == 1);
    mdac_env* cur = nullptr;
    OK(mdac_env_current(&cur));
    REQUIRE(cur == envs.env);
    require_same(coeffs(y), coeffs(x));
    H z;
    OK(mdac_nda_import(envs.env, y, z.out()));
    REQUIRE(mdac_nda_env(z) == envs.env);
    require_same(coeffs(z), coeffs(x));

    C c, cb;
    OK(mdac_cnda_new(x, z, c.out()));
    OK(mdac_cnda_import(b, c, cb.out()));
    REQUIRE(mdac_cnda_env(cb) == b);
    H re, im;
    OK(mdac_cnda_real(cb, re.out()));
    OK(mdac_cnda_imag(cb, im.out()));
    require_same(coeffs(re), coeffs(x));
    require_same(coeffs(im), coeffs(x));

    H bad;
    REQUIRE(mdac_nda_import(small, x, bad.out()) == MDAC_ERR_VALUE);   // different layout
    C cbad;
    REQUIRE(mdac_cnda_import(small, c, cbad.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_nda_import(nullptr, x, bad.out()) == MDAC_ERR_ENV);
    mdac_nda_free(y.p);
    y.p = nullptr;
    mdac_cnda_free(cb.p);
    cb.p = nullptr;
    OK(mdac_env_close(b));
    REQUIRE(mdac_nda_import(b, x, bad.out()) == MDAC_ERR_ENV);         // closed target
    REQUIRE(mdac_cnda_import(b, c, cbad.out()) == MDAC_ERR_ENV);
    OK(mdac_env_close(small));
}

TEST_CASE("capi: importing from a closed env fails", "[capi][multienv]") {
    Envs envs;
    mdac_env* b = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 20, 0, &b));
    H x;
    OK(mdac_nda_var(b, 0, x.out()));
    OK(mdac_env_close(b));
    H y;
    REQUIRE(mdac_nda_import(envs.env, x, y.out()) == MDAC_ERR_ENV);
    REQUIRE(mdac_nda_env(x) == b);
}

// ---- Symbolic: Expr and SDA (T5.1, T5.2) ---------------------------------------

TEST_CASE("capi: mdac_has_symbolic", "[capi][sym]") {
#ifdef DA_WITH_SYMBOLIC
    REQUIRE(mdac_has_symbolic() == 1);
#else
    REQUIRE(mdac_has_symbolic() == 0);
    mdac_expr* x = nullptr;
    REQUIRE(mdac_expr_new_d(1.0, &x) == MDAC_ERR_UNSUPPORTED);
    REQUIRE(std::string(mdac_last_error()).find("symbolic") != std::string::npos);
    Envs envs;
    mdac_csda* c = nullptr;
    REQUIRE(mdac_csda_new_z(envs.env, 1.0, 0.0, &c) == MDAC_ERR_UNSUPPORTED);
    mdac_csdalist* l = nullptr;
    REQUIRE(mdac_csdalist_new(&l) == MDAC_ERR_UNSUPPORTED);
    REQUIRE(mdac_csdalist_length(l) == 0);
    REQUIRE(mdac_csda_env(c) == nullptr);
    mdac_csda_free(c);
    mdac_sda* sp = nullptr;
    REQUIRE(mdac_nda_promote_to(envs.env, nullptr, &sp) == MDAC_ERR_UNSUPPORTED);
    REQUIRE(mdac_sda_import(envs.env, nullptr, &sp) == MDAC_ERR_UNSUPPORTED);
    REQUIRE(mdac_csda_import(envs.env, nullptr, &c) == MDAC_ERR_UNSUPPORTED);
#endif
}

#ifdef DA_WITH_SYMBOLIC

#include <symengine/parser.h>

namespace {

using SymEngine::Expression;

// Owns a C API Expr.
struct EX {
    mdac_expr* p = nullptr;
    EX() = default;
    explicit EX(mdac_expr* q) : p(q) {}
    EX(const EX&) = delete;
    EX& operator=(const EX&) = delete;
    ~EX() { mdac_expr_free(p); }
    operator mdac_expr*() const { return p; }
    mdac_expr** out() { return &p; }
};

// Owns a C API SDA.
struct S {
    mdac_sda* p = nullptr;
    S() = default;
    explicit S(mdac_sda* q) : p(q) {}
    S(const S&) = delete;
    S& operator=(const S&) = delete;
    ~S() { mdac_sda_free(p); }
    operator mdac_sda*() const { return p; }
    mdac_sda** out() { return &p; }
};

// Owns a C API SDA list.
struct SL {
    mdac_sdalist* p = nullptr;
    SL() = default;
    SL(const SL&) = delete;
    SL& operator=(const SL&) = delete;
    ~SL() { mdac_sdalist_free(p); }
    operator mdac_sdalist*() const { return p; }
    mdac_sdalist** out() { return &p; }
};

std::string str(const mdac_expr* x) {
    size_t n = 0;
    OK(mdac_expr_to_string(x, nullptr, 0, &n));
    std::string s(n, '\0');
    OK(mdac_expr_to_string(x, s.data(), n, &n));
    s.resize(n - 1);
    return s;
}

std::string str(const Expression& e) {
    std::ostringstream os;
    os << e;
    return os.str();
}

mdac_expr* parse(const char* s) {
    mdac_expr* x = nullptr;
    OK(mdac_expr_parse(s, &x));
    return x;
}

std::vector<std::string> scoeffs(const mdac_sda* v) {
    size_t n = 0;
    OK(mdac_sda_coeffs(v, nullptr, 0, &n));
    std::vector<mdac_expr*> buf(n);
    OK(mdac_sda_coeffs(v, buf.data(), n, &n));
    std::vector<std::string> out;
    for (mdac_expr* x : buf) {
        out.push_back(str(x));
        mdac_expr_free(x);
    }
    return out;
}

std::vector<std::string> scoeffs(const da::SDA& v) {
    const Expression* p = v.env_->pool<Expression>().slot(v.slot_);
    std::vector<std::string> out;
    for (unsigned i = 0; i < v.length(); ++i) out.push_back(str(p[i]));
    return out;
}

void require_same(const mdac_sda* c, const da::SDA& r) { REQUIRE(scoeffs(c) == scoeffs(r)); }

// The test inputs: coefficients (as text) of the monomials in exps, in both APIs.
const std::vector<std::vector<int>> kExps = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 2}};

mdac_sda* sfrom(mdac_env* e, const std::vector<const char*>& cs) {
    mdac_sda* v = nullptr;
    OK(mdac_sda_new(e, nullptr, &v));
    for (size_t i = 0; i < cs.size(); ++i) {
        EX x(parse(cs[i]));
        OK(mdac_sda_set_coeff(v, kExps[i].data(), kExps[i].size(), x));
    }
    return v;
}

da::SDA sref(const std::vector<const char*>& cs) {
    da::SDA v;
    for (size_t i = 0; i < cs.size(); ++i) v.set_element(kExps[i], Expression(SymEngine::parse(cs[i])));
    return v;
}

template <class F, class G>
void sboth(mdac_sda* o, F alloc, G into, const da::SDA& want) {
    {
        S r;
        OK(alloc(r.out()));
        require_same(r, want);
    }
    OK(into(o));
    require_same(o, want);
}

#define SALLOC(expr) [&](mdac_sda** r_) { return expr; }
#define SINTO(expr) [&](mdac_sda* t_) { return expr; }

#define SDA_OP(name, op)                                                                          \
    SECTION(#name) {                                                                              \
        sboth(o, SALLOC(mdac_sda_##name(a, b, r_)), SINTO(mdac_sda_##name##_into(t_, a, b)), ra op rb); \
        sboth(o, SALLOC(mdac_sda_##name##_n(a, x, r_)), SINTO(mdac_sda_##name##_n_into(t_, a, x)), ra op rx); \
        sboth(o, SALLOC(mdac_sda_n##name(x, a, r_)), SINTO(mdac_sda_n##name##_into(t_, x, a)), rx op ra); \
        sboth(o, SALLOC(mdac_sda_##name##_e(a, ex, r_)), SINTO(mdac_sda_##name##_e_into(t_, a, ex)), ra op rex); \
        sboth(o, SALLOC(mdac_sda_e##name(ex, a, r_)), SINTO(mdac_sda_e##name##_into(t_, ex, a)), rex op ra); \
        sboth(o, SALLOC(mdac_sda_##name##_d(a, d, r_)), SINTO(mdac_sda_##name##_d_into(t_, a, d)), ra op d); \
        sboth(o, SALLOC(mdac_sda_d##name(d, a, r_)), SINTO(mdac_sda_d##name##_into(t_, d, a)), d op ra); \
        sboth(o, SALLOC(mdac_nda_##name##_e(x, ex, r_)), SINTO(mdac_nda_##name##_e_into(t_, x, ex)), rx op rex); \
        sboth(o, SALLOC(mdac_nda_e##name(ex, x, r_)), SINTO(mdac_nda_e##name##_into(t_, ex, x)), rex op rx); \
        S a1(sfrom(e, ca)), b1(sfrom(e, cb));                                                     \
        OK(mdac_sda_##name##_into(a1, a1, b));                                                    \
        require_same(a1, ra op rb);                                                               \
        OK(mdac_sda_##name##_into(b1, a, b1));                                                    \
        require_same(b1, ra op rb);                                                               \
    }

const std::vector<const char*> ca = {"1 + a", "b", "2", "a*b", "1/3"};
const std::vector<const char*> cb = {"2", "a", "0", "0", "b**2"};

} // namespace

TEST_CASE("capi: Expr lifecycle and conversions", "[capi][sym][expr]") {
    EX d, i, p, s, c, big;
    OK(mdac_expr_new_d(0.25, d.out()));
    OK(mdac_expr_new_i(3, i.out()));
    OK(mdac_expr_new_i(INT64_MIN, big.out()));
    OK(mdac_expr_symbol("a", s.out()));
    p.p = parse("a + 2*b");
    OK(mdac_expr_copy(p, c.out()));
    REQUIRE(str(d) == str(Expression(0.25)));
    REQUIRE(str(i) == "3");
    REQUIRE(str(big) == "-9223372036854775808");
    REQUIRE(str(s) == "a");
    REQUIRE(str(c) == str(Expression(SymEngine::parse("a + 2*b"))));

    size_t n = 0;
    char buf[4];
    OK(mdac_expr_to_string(p, buf, sizeof buf, &n));
    REQUIRE(n == str(p).size() + 1);
    REQUIRE(std::string(buf) == str(p).substr(0, 3));

    double v = 0;
    OK(mdac_expr_to_double(d, &v));
    REQUIRE(v == 0.25);
    EX sin1(parse("sin(1) + 2"));
    OK(mdac_expr_to_double(sin1, &v));
    REQUIRE(v == Approx(std::sin(1.0) + 2));
    REQUIRE(mdac_expr_to_double(p, &v) == MDAC_ERR_VALUE);
    REQUIRE(std::string(mdac_last_error()).size() > 0);

    EX bad;
    REQUIRE(mdac_expr_parse("a +* b", bad.out()) == MDAC_ERR_VALUE);
    REQUIRE(bad.p == nullptr);
}

TEST_CASE("capi: Expr arithmetic and methods equal C++", "[capi][sym][expr]") {
    const Expression ra(SymEngine::symbol("a")), rb(SymEngine::symbol("b"));
    EX a, b;
    OK(mdac_expr_symbol("a", a.out()));
    OK(mdac_expr_symbol("b", b.out()));
    const double d = 1.5;
    auto check = [](mdac_status st, EX& r, const Expression& want) {
        REQUIRE(st == MDAC_OK);
        REQUIRE(str(r) == str(want));
    };
#define EXPR_OP(name, expr_ab, expr_ad, expr_da)                     \
    {                                                                \
        EX r1, r2, r3;                                                \
        check(mdac_expr_##name(a, b, r1.out()), r1, expr_ab);        \
        check(mdac_expr_##name##_d(a, d, r2.out()), r2, expr_ad);    \
        check(mdac_expr_d##name(d, a, r3.out()), r3, expr_da);       \
    }
    EXPR_OP(add, ra + rb, ra + Expression(d), Expression(d) + ra)
    EXPR_OP(sub, ra - rb, ra - Expression(d), Expression(d) - ra)
    EXPR_OP(mul, ra * rb, ra * Expression(d), Expression(d) * ra)
    EXPR_OP(div, ra / rb, ra / Expression(d), Expression(d) / ra)
    EXPR_OP(pow, SymEngine::pow(ra, rb), SymEngine::pow(ra, Expression(d)),
            SymEngine::pow(Expression(d), ra))
#undef EXPR_OP
    EX ng;
    check(mdac_expr_neg(a, ng.out()), ng, -ra);

    EX s1, s2;
    OK(mdac_expr_add(a, b, s1.out()));
    OK(mdac_expr_add(b, a, s2.out()));
    int eq = 0;
    OK(mdac_expr_eq(s1, s2, &eq));
    REQUIRE(eq == 1);
    OK(mdac_expr_eq(a, b, &eq));
    REQUIRE(eq == 0);
    uint64_t h1 = 0, h2 = 0;
    OK(mdac_expr_hash(s1, &h1));
    OK(mdac_expr_hash(s2, &h2));
    REQUIRE(h1 == h2);

    EX sq(parse("(a + b)**2")), ex, df, sm, sb;
    check(mdac_expr_expand(sq, ex.out()), ex, SymEngine::expand(Expression(SymEngine::parse("(a + b)**2"))));
    check(mdac_expr_diff(sq, a, df.out()), df, Expression(SymEngine::parse("2*(a + b)")));
    EX notsym;
    REQUIRE(mdac_expr_diff(sq, s1, notsym.out()) == MDAC_ERR_VALUE);
    EX twice(parse("a*2/2"));
    check(mdac_expr_simplify(twice, sm.out()), sm, ra);
    EX one;
    OK(mdac_expr_new_d(1.0, one.out()));
    const mdac_expr* keys[] = {a, b};
    const mdac_expr* vals[] = {one, a};
    check(mdac_expr_subs(sq, 2, keys, vals, sb.out()), sb, Expression(SymEngine::parse("(1.0 + a)**2")));

    int z = -1;
    EX diff0;
    OK(mdac_expr_sub(ex, sq, diff0.out()));
    OK(mdac_expr_is_zero(diff0, &z));
    REQUIRE(z == 1);
    OK(mdac_expr_is_zero(sq, &z));
    REQUIRE(z == 0);

    size_t n = 0;
    OK(mdac_expr_free_symbols(sq, nullptr, 0, &n));
    REQUIRE(n == 2);
    mdac_expr* syms[2] = {nullptr, nullptr};
    OK(mdac_expr_free_symbols(sq, syms, 1, &n));
    REQUIRE(n == 2);
    REQUIRE(syms[1] == nullptr);
    mdac_expr_free(syms[0]);
    OK(mdac_expr_free_symbols(sq, syms, 2, &n));
    std::set<std::string> names = {str(syms[0]), str(syms[1])};
    REQUIRE(names == std::set<std::string>{"a", "b"});
    mdac_expr_free(syms[0]);
    mdac_expr_free(syms[1]);
}

TEST_CASE("capi: SDA lifecycle and inspection", "[capi][sym][sda]") {
    Envs envs;
    mdac_env* e = envs.env;
    S zero, c, cd, v, p, cp;
    OK(mdac_sda_new(e, nullptr, zero.out()));
    REQUIRE(scoeffs(zero) == std::vector<std::string>{"0"});
    EX ex(parse("a"));
    OK(mdac_sda_new(e, ex, c.out()));
    REQUIRE(scoeffs(c) == std::vector<std::string>{"a"});
    OK(mdac_sda_new_d(e, 1.5, cd.out()));
    require_same(cd, da::SDA(1.5));
    OK(mdac_sda_var(e, 1, v.out()));
    require_same(v, da::promote(da::da_base(1)));
    S bad;
    REQUIRE(mdac_sda_var(e, kNvars, bad.out()) == MDAC_ERR_INDEX);

    const da::NDA rx = nda(1.5, 2.0);
    H x(from(e, coeffs(rx)));
    OK(mdac_sda_promote(x, p.out()));
    require_same(p, da::promote(rx));

    S a(sfrom(e, ca));
    const da::SDA ra = sref(ca);
    require_same(a, ra);
    OK(mdac_sda_copy(a, cp.out()));
    REQUIRE(cp.p != a.p);
    require_same(cp, ra);
    REQUIRE(mdac_sda_env(a) == e);

    EX con;
    OK(mdac_sda_con(a, con.out()));
    REQUIRE(str(con) == "1 + a");
    size_t len = 0, nt = 0;
    OK(mdac_sda_length(a, &len));
    OK(mdac_sda_nterms(a, &nt));
    REQUIRE(len == ra.length());
    REQUIRE(nt == static_cast<size_t>(const_cast<da::SDA&>(ra).n_element()));

    EX k;
    const int ab[] = {1, 1};
    OK(mdac_sda_coeff(a, ab, 2, k.out()));
    REQUIRE(str(k) == "a*b");
    const int neg[] = {-1, 0, 0};
    EX kneg;
    REQUIRE(mdac_sda_coeff(a, neg, 3, kneg.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_sda_set_coeff(a, neg, 3, ex) == MDAC_ERR_VALUE);

    int exps[kNvars] = {};
    EX t;
    OK(mdac_sda_index_term(a, 1, exps, t.out()));
    std::vector<unsigned> rexps(kNvars);
    Expression rt;
    da::detail::ad_elem(ra.env_->layout(), ra.env_->pool<Expression>(), ra.slot_, 2, rexps.data(), rt);
    REQUIRE(str(t) == str(rt));
    REQUIRE(std::vector<int>(exps, exps + kNvars) == std::vector<int>(rexps.begin(), rexps.end()));
    EX tbad;
    REQUIRE(mdac_sda_index_term(a, da::NDA::full_length(), exps, tbad.out()) == MDAC_ERR_INDEX);

    // The operator<< text; its first line names the slot.
    size_t n = 0;
    OK(mdac_sda_to_string(a, nullptr, 0, &n));
    std::string text(n, '\0');
    OK(mdac_sda_to_string(a, text.data(), n, &n));
    text.resize(n - 1);
    std::ostringstream os;
    os << ra;
    REQUIRE(text.substr(text.find('\n')) == os.str().substr(os.str().find('\n')));

    int z = -1;
    OK(mdac_sda_iszero(a, &z));
    REQUIRE(z == 0);
    S diff;
    OK(mdac_sda_sub(a, a, diff.out()));
    OK(mdac_sda_iszero(diff, &z));
    REQUIRE(z == 1);
    S u(sfrom(e, {"1", "a - a"}));
    OK(mdac_sda_clean(u));
    OK(mdac_sda_length(u, &len));
    REQUIRE(len == 1);
    OK(mdac_sda_set_con(cp, ex));
    REQUIRE(scoeffs(cp) == std::vector<std::string>{"a"});
    OK(mdac_sda_reset(cp));
    OK(mdac_sda_iszero(cp, &z));
    REQUIRE(z == 1);
}

TEST_CASE("capi: SDA arithmetic equals C++", "[capi][sym][sda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::SDA ra = sref(ca), rb = sref(cb);
    const da::NDA rx = nda(0.9, -1.0);
    const Expression rex = SymEngine::parse("b + 1/2");
    const double d = 1.25;
    S a(sfrom(e, ca)), b(sfrom(e, cb)), o(sfrom(e, {"7"}));
    H x(from(e, coeffs(rx)));
    EX ex(parse("b + 1/2"));

    SDA_OP(add, +)
    SDA_OP(sub, -)
    SDA_OP(mul, *)
    SDA_OP(div, /)
    SECTION("neg and pow") {
        sboth(o, SALLOC(mdac_sda_neg(a, r_)), SINTO(mdac_sda_neg_into(t_, a)), -ra);
        sboth(o, SALLOC(mdac_sda_pow_i(a, 3, r_)), SINTO(mdac_sda_pow_i_into(t_, a, 3)), da::pow(ra, 3));
        sboth(o, SALLOC(mdac_sda_pow_d(a, 0.3, r_)), SINTO(mdac_sda_pow_d_into(t_, a, 0.3)),
              da::pow(ra, 0.3));
        OK(mdac_sda_neg_into(a, a));
        require_same(a, -ra);
    }
    SECTION("division by zero") {
        S r;
        REQUIRE(mdac_sda_div_d(a, 0.0, r.out()) == MDAC_ERR_VALUE);
        REQUIRE(mdac_sda_div_d_into(o, a, 0.0) == MDAC_ERR_VALUE);
    }
}

TEST_CASE("capi: SDA math functions equal C++", "[capi][sym][sda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const da::SDA ra = sref({"1/2 + a", "1", "b"});
    S a(sfrom(e, {"1/2 + a", "1", "b"})), o(sfrom(e, {"7"}));
#define SDA_FN(f) sboth(o, SALLOC(mdac_sda_##f(a, r_)), SINTO(mdac_sda_##f##_into(t_, a)), da::f(ra));
    SDA_FN(sqrt) SDA_FN(exp) SDA_FN(log) SDA_FN(sin) SDA_FN(cos) SDA_FN(tan) SDA_FN(asin)
    SDA_FN(acos) SDA_FN(atan) SDA_FN(sinh) SDA_FN(cosh) SDA_FN(tanh) SDA_FN(erf)
#undef SDA_FN
}

TEST_CASE("capi: SDA simplify, expand, subs and evaluate equal C++", "[capi][sym][sda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const std::vector<const char*> cs = {"(a + 1)**2", "a*2/2 + b - b", "sin(a)**2"};
    const da::SDA rs = sref(cs);
    S s(sfrom(e, cs)), sp, ep, sb;
    OK(mdac_sda_simplify(s, sp.out()));
    da::SDA want = rs;
    want.simplify();
    require_same(sp, want);
    OK(mdac_sda_expand(s, ep.out()));
    REQUIRE(scoeffs(ep)[0] == str(SymEngine::expand(Expression(SymEngine::parse("(a + 1)**2")))));
    require_same(s, rs);                                       // unchanged

    EX a(parse("a")), b(parse("b")), two;
    OK(mdac_expr_new_i(2, two.out()));
    const mdac_expr* keys[] = {a};
    const mdac_expr* vals[] = {two};
    OK(mdac_sda_subs(s, 1, keys, vals, sb.out()));
    REQUIRE(scoeffs(sb)[0] == "9");

    const mdac_expr* ks[] = {a, b};
    const double xs[] = {0.3, -1.25};
    H r;
    OK(mdac_sda_evaluate(s, 2, ks, xs, r.out()));
    SymEngine::vec_basic syms = {SymEngine::symbol("a"), SymEngine::symbol("b")};
    require_same(coeffs(r), coeffs(da::evaluate(rs, syms, {0.3, -1.25})));
    H r2;
    REQUIRE(mdac_sda_evaluate(s, 1, ks + 1, xs, r2.out()) == MDAC_ERR_VALUE);   // a has no value
}

TEST_CASE("capi: SDA lists and algorithms equal C++", "[capi][sym][sda][algo]") {
    Envs envs;
    mdac_env* e = envs.env;
    da::SDA ra = sref(ca), rb = sref(cb);
    S a(sfrom(e, ca)), b(sfrom(e, cb));

    SECTION("lists") {
        SL l;
        OK(mdac_sdalist_new(l.out()));
        REQUIRE(mdac_sdalist_length(l) == 0);
        REQUIRE(mdac_sdalist_env(l) == nullptr);
        OK(mdac_sdalist_push(l, a));
        OK(mdac_sdalist_push(l, b));
        REQUIRE(mdac_sdalist_length(l) == 2);
        REQUIRE(mdac_sdalist_env(l) == e);
        S g;
        OK(mdac_sdalist_get(l, 1, g.out()));
        require_same(g, rb);
        OK(mdac_sdalist_set(l, 1, a));
        S g2;
        OK(mdac_sdalist_get(l, 1, g2.out()));
        require_same(g2, ra);
        S g3;
        REQUIRE(mdac_sdalist_get(l, 2, g3.out()) == MDAC_ERR_INDEX);
        const mdac_sda* vs[] = {a, b};
        SL f;
        OK(mdac_sdalist_from(vs, 2, f.out()));
        REQUIRE(mdac_sdalist_length(f) == 2);
    }
    SECTION("der, integ and substitute") {
        S r1, r2, r3, r4, r5, bad;
        OK(mdac_sda_der(a, 0, r1.out()));
        require_same(r1, da::da_der(ra, 0));
        OK(mdac_sda_integ(a, 1, r2.out()));
        require_same(r2, da::da_int(ra, 1));
        REQUIRE(mdac_sda_der(a, kNvars, bad.out()) == MDAC_ERR_INDEX);
        da::SDA w;
        OK(mdac_sda_substitute_d(a, 0, 2.0, r3.out()));
        da::da_substitute_const(ra, 0, 2.0, w);
        require_same(r3, w);
        OK(mdac_sda_substitute(a, 1, b, r4.out()));
        da::da_substitute(ra, 1, rb, w);
        require_same(r4, w);
        const mdac_sda* vs[] = {b, a};
        SL l;
        OK(mdac_sdalist_from(vs, 2, l.out()));
        const unsigned ids[] = {0, 2};
        OK(mdac_sda_substitute_multi(a, ids, 2, l, r5.out()));
        std::vector<unsigned> rids = {0, 2};
        std::vector<da::SDA> rl = {rb, ra};
        da::da_substitute(ra, rids, rl, w);
        require_same(r5, w);
        const unsigned dup[] = {1, 1};
        REQUIRE(mdac_sda_substitute_multi(a, dup, 2, l, bad.out()) == MDAC_ERR_VALUE);

        const mdac_sda* ms[] = {a, b};
        SL m, o;
        OK(mdac_sdalist_from(ms, 2, m.out()));
        OK(mdac_sdalist_substitute(m, ids, 2, l, o.out()));
        std::vector<da::SDA> rm = {ra, rb}, ro(2);
        da::da_substitute(rm, rids, rl, ro);
        for (size_t i = 0; i < 2; ++i) {
            S g;
            OK(mdac_sdalist_get(o, i, g.out()));
            require_same(g, ro[i]);
        }
    }
    SECTION("compose") {
        const da::SDA rc = sref({"1/2", "a"});
        S c(sfrom(e, {"1/2", "a"}));
        const mdac_sda* ms[] = {a, b};
        const mdac_sda* vs[] = {c, a, b};
        SL m, v, o, v2;
        OK(mdac_sdalist_from(ms, 2, m.out()));
        OK(mdac_sdalist_from(vs, 3, v.out()));
        OK(mdac_sdalist_compose(m, v, o.out()));
        std::vector<da::SDA> rm = {ra, rb}, rv = {rc, ra, rb}, ro(2);
        da::da_composition(rm, rv, ro);
        for (size_t i = 0; i < 2; ++i) {
            S g;
            OK(mdac_sdalist_get(o, i, g.out()));
            require_same(g, ro[i]);
        }
        OK(mdac_sdalist_from(vs, 2, v2.out()));
        SL bad;
        REQUIRE(mdac_sdalist_compose(m, v2, bad.out()) == MDAC_ERR_VALUE);

        const double pt[] = {0.5, -1.0, 2.0};
        mdac_expr* vals[2] = {nullptr, nullptr};
        OK(mdac_sdalist_compose_d(m, pt, 3, vals));
        std::vector<double> rpt(pt, pt + 3);
        std::vector<Expression> rvals;
        da::da_composition(rm, rpt, rvals);
        for (size_t i = 0; i < 2; ++i) {
            REQUIRE(str(vals[i]) == str(rvals[i]));
            mdac_expr_free(vals[i]);
        }
        REQUIRE(mdac_sdalist_compose_d(m, pt, 2, vals) == MDAC_ERR_VALUE);
    }
}

TEST_CASE("capi: SDA envs: mixing fails, freeing after mdac_clear is safe", "[capi][sym][sda]") {
    mdac_sda* s = nullptr;
    {
        Envs envs;
        OK(mdac_sda_var(envs.env, 0, &s));
        mdac_env* other = nullptr;
        OK(mdac_env_make(kOrder, kNvars, 50, 0, &other));
        S t;
        OK(mdac_sda_new_d(other, 1.0, t.out()));
        H x;
        OK(mdac_nda_new(other, 2.0, x.out()));
        S r;
        REQUIRE(mdac_sda_add(s, t, r.out()) == MDAC_ERR_ENV);
        REQUIRE(mdac_sda_mul_n(s, x, r.out()) == MDAC_ERR_ENV);
        REQUIRE(mdac_sda_add_into(t, s, s) == MDAC_ERR_ENV);
        mdac_sda_free(t.p);
        t.p = nullptr;
        mdac_nda_free(x.p);
        x.p = nullptr;
        OK(mdac_env_close(other));
    }
    S r;
    REQUIRE(mdac_sda_exp(s, r.out()) == MDAC_ERR_ENV);
    mdac_sda_free(s);
}

// ---- T6.1 ------------------------------------------------------------------

namespace {

// Owns a C API CSDA.
struct CS {
    mdac_csda* p = nullptr;
    CS() = default;
    explicit CS(mdac_csda* q) : p(q) {}
    CS(const CS&) = delete;
    CS& operator=(const CS&) = delete;
    ~CS() { mdac_csda_free(p); }
    operator mdac_csda*() const { return p; }
    mdac_csda** out() { return &p; }
};

// Owns a C API CSDA list.
struct CSL {
    mdac_csdalist* p = nullptr;
    CSL() = default;
    CSL(const CSL&) = delete;
    CSL& operator=(const CSL&) = delete;
    ~CSL() { mdac_csdalist_free(p); }
    operator mdac_csdalist*() const { return p; }
    mdac_csdalist** out() { return &p; }
};

using RCSDA = std::complex<da::SDA>;

mdac_csda* csfrom(mdac_env* e, const std::vector<const char*>& re, const std::vector<const char*>& im) {
    S r(sfrom(e, re)), i(sfrom(e, im));
    mdac_csda* c = nullptr;
    OK(mdac_csda_new(r, i, &c));
    return c;
}

RCSDA csref(const std::vector<const char*>& re, const std::vector<const char*>& im) {
    return RCSDA(sref(re), sref(im));
}

void require_same(const mdac_csda* c, const RCSDA& r) {
    S re, im;
    OK(mdac_csda_real(c, re.out()));
    OK(mdac_csda_imag(c, im.out()));
    require_same(re, get_real(r));
    require_same(im, get_imag(r));
}

void require_same(const mdac_cnda* c, const da::CNDA& r, double tol) {
    H re, im;
    OK(mdac_cnda_real(c, re.out()));
    OK(mdac_cnda_imag(c, im.out()));
    std::vector<double> a = coeffs(re), b = coeffs(get_real(r)), x = coeffs(im), y = coeffs(get_imag(r));
    a.resize(std::max(a.size(), b.size()));
    b.resize(a.size());
    x.resize(std::max(x.size(), y.size()));
    y.resize(x.size());
    for (size_t i = 0; i < a.size(); ++i) REQUIRE(a[i] == Approx(b[i]).margin(tol));
    for (size_t i = 0; i < x.size(); ++i) REQUIRE(x[i] == Approx(y[i]).margin(tol));
}

// The allocating form and the _into form (into o, whose slots it keeps) both give want.
template <class F, class G>
void csboth(mdac_env* e, mdac_csda* o, F alloc, G into, const RCSDA& want) {
    {
        CS r;
        OK(alloc(r.out()));
        require_same(r, want);
    }
    const unsigned n = env_count(e);
    OK(into(o));
    REQUIRE(env_count(e) == n);
    require_same(o, want);
}

#define CSALLOC(expr) [&](mdac_csda** r_) { return expr; }
#define CSINTO(expr) [&](mdac_csda* t_) { return expr; }

#define CSDA_OP(name, op)                                                                         \
    SECTION(#name) {                                                                              \
        csboth(e, o, CSALLOC(mdac_csda_##name(a, b, r_)), CSINTO(mdac_csda_##name##_into(t_, a, b)), ra op rb); \
        csboth(e, o, CSALLOC(mdac_csda_##name##_s(a, s, r_)), CSINTO(mdac_csda_##name##_s_into(t_, a, s)), ra op rs); \
        csboth(e, o, CSALLOC(mdac_csda_s##name(s, a, r_)), CSINTO(mdac_csda_s##name##_into(t_, s, a)), rs op ra); \
        csboth(e, o, CSALLOC(mdac_csda_##name##_e(a, ex, r_)), CSINTO(mdac_csda_##name##_e_into(t_, a, ex)), ra op rex); \
        csboth(e, o, CSALLOC(mdac_csda_e##name(ex, a, r_)), CSINTO(mdac_csda_e##name##_into(t_, ex, a)), rex op ra); \
        csboth(e, o, CSALLOC(mdac_csda_##name##_d(a, d, r_)), CSINTO(mdac_csda_##name##_d_into(t_, a, d)), ra op d); \
        csboth(e, o, CSALLOC(mdac_csda_d##name(d, a, r_)), CSINTO(mdac_csda_d##name##_into(t_, d, a)), d op ra); \
        csboth(e, o, CSALLOC(mdac_csda_##name##_z(a, z.real(), z.imag(), r_)),                      \
               CSINTO(mdac_csda_##name##_z_into(t_, a, z.real(), z.imag())), ra op z);           \
        csboth(e, o, CSALLOC(mdac_csda_z##name(z.real(), z.imag(), a, r_)),                         \
               CSINTO(mdac_csda_z##name##_into(t_, z.real(), z.imag(), a)), z op ra);            \
        csboth(e, o, CSALLOC(mdac_sda_##name##_z(s, z.real(), z.imag(), r_)),                       \
               CSINTO(mdac_sda_##name##_z_into(t_, s, z.real(), z.imag())), rs op z);            \
        csboth(e, o, CSALLOC(mdac_sda_z##name(z.real(), z.imag(), s, r_)),                          \
               CSINTO(mdac_sda_z##name##_into(t_, z.real(), z.imag(), s)), z op rs);             \
        CS a1(csfrom(e, ca, cb)), b1(csfrom(e, cc, cd));                                          \
        OK(mdac_csda_##name##_into(a1, a1, b));                                                   \
        require_same(a1, ra op rb);                                                               \
        OK(mdac_csda_##name##_into(b1, a, b1));                                                   \
        require_same(b1, ra op rb);                                                               \
    }

const std::vector<const char*> cc = {"1/2 + b", "1", "a"};
const std::vector<const char*> cd = {"1/4", "0", "2*b"};

// The operator<< text of a CSDA without the line after "... part:", which
// names the part's slot.
std::string without_slot_lines(const std::string& text) {
    std::istringstream is(text);
    std::string line, out;
    bool skip = false;
    while (std::getline(is, line)) {
        if (!skip) out += line + "\n";
        skip = line.size() > 5 && line.compare(line.size() - 5, 5, "part:") == 0;
    }
    return out;
}

// Every element of the C list equals the C++ one.
template <class L>
void require_same(const mdac_csdalist* l, const L& r) {
    REQUIRE(mdac_csdalist_length(l) == r.size());
    for (size_t i = 0; i < r.size(); ++i) {
        CS g;
        OK(mdac_csdalist_get(l, i, g.out()));
        require_same(g, r[i]);
    }
}

} // namespace

TEST_CASE("capi: CSDA lifecycle and parts", "[capi][sym][csda]") {
    Envs envs;
    mdac_env* e = envs.env;
    const RCSDA ra = csref(ca, cb);
    CS a(csfrom(e, ca, cb)), re_only, z, p, cp;
    require_same(a, ra);
    S sa(sfrom(e, ca));
    OK(mdac_csda_new(sa, nullptr, re_only.out()));
    require_same(re_only, RCSDA(sref(ca), da::SDA()));
    OK(mdac_csda_new_z(e, 1.5, -2.0, z.out()));
    require_same(z, RCSDA(da::SDA(1.5), da::SDA(-2.0)));

    const da::CNDA rn(nda(1.5, 2.0), nda(-0.5, 1.0));
    H nr(from(e, coeffs(get_real(rn)))), ni(from(e, coeffs(get_imag(rn))));
    C n;
    OK(mdac_cnda_new(nr, ni, n.out()));
    OK(mdac_csda_promote(n, p.out()));
    require_same(p, da::promote(rn));

    OK(mdac_csda_copy(a, cp.out()));
    REQUIRE(cp.p != a.p);
    require_same(cp, ra);
    REQUIRE(mdac_csda_env(a) == e);
    S sb(sfrom(e, cb));
    OK(mdac_csda_set_real(cp, sb));
    OK(mdac_csda_set_imag(cp, sa));
    require_same(cp, RCSDA(sref(cb), sref(ca)));
    require_same(a, ra);                                       // the copy was independent

    // The operator<< text; the first line of each part's table names its slot.
    size_t len = 0;
    OK(mdac_csda_to_string(a, nullptr, 0, &len));
    std::string text(len, '\0');
    OK(mdac_csda_to_string(a, text.data(), len, &len));
    text.resize(len - 1);
    std::ostringstream os;
    os << ra;
    REQUIRE(text.find("Imaginary part") != std::string::npos);
    REQUIRE(without_slot_lines(text) == without_slot_lines(os.str()));
    char small[5];
    OK(mdac_csda_to_string(a, small, sizeof small, &len));
    REQUIRE(std::string(small) == text.substr(0, 4));
}

TEST_CASE("capi: CSDA arithmetic equals C++", "[capi][sym][csda]") {
    Envs envs(400);
    mdac_env* e = envs.env;
    const RCSDA ra = csref(ca, cb), rb = csref(cc, cd);
    const da::SDA rs = sref({"1/3 + a", "b"});
    const Expression rex = SymEngine::parse("b + 1/2");
    const double d = 1.25;
    const std::complex<double> z(0.5, -2.0);
    CS a(csfrom(e, ca, cb)), b(csfrom(e, cc, cd)), o(csfrom(e, {"7"}, {"1"}));
    S s(sfrom(e, {"1/3 + a", "b"}));
    EX ex(parse("b + 1/2"));

    CSDA_OP(add, +)
    CSDA_OP(sub, -)
    CSDA_OP(mul, *)
    CSDA_OP(div, /)
    SECTION("neg") {
        csboth(e, o, CSALLOC(mdac_csda_neg(a, r_)), CSINTO(mdac_csda_neg_into(t_, a)), -ra);
        OK(mdac_csda_neg_into(a, a));
        require_same(a, -ra);
    }
    SECTION("division by zero") {
        CS r;
        REQUIRE(mdac_csda_div_d(a, 0.0, r.out()) == MDAC_ERR_VALUE);
        REQUIRE(mdac_csda_div_d_into(o, a, 0.0) == MDAC_ERR_VALUE);
    }
}

TEST_CASE("capi: CSDA math functions, abs and evaluate equal C++", "[capi][sym][csda]") {
    Envs envs(400);
    mdac_env* e = envs.env;
    // Numeric constant parts keep the symbolic series small.
    const std::vector<const char*> re = {"1/2", "a"}, im = {"1/3", "0", "b"};
    const RCSDA ra = csref(re, im);
    CS a(csfrom(e, re, im)), o(csfrom(e, {"7"}, {"1"}));
#define CSDA_FN(f) csboth(e, o, CSALLOC(mdac_csda_##f(a, r_)), CSINTO(mdac_csda_##f##_into(t_, a)), da::f(ra));
    CSDA_FN(sqrt) CSDA_FN(exp) CSDA_FN(log) CSDA_FN(asin) CSDA_FN(acos) CSDA_FN(atan)
    CSDA_FN(asinh) CSDA_FN(acosh) CSDA_FN(atanh)
#undef CSDA_FN
    csboth(e, o, CSALLOC(mdac_csda_pow_i(a, 3, r_)), CSINTO(mdac_csda_pow_i_into(t_, a, 3)),
           da::pow(ra, 3));
    csboth(e, o, CSALLOC(mdac_csda_pow_d(a, 0.5, r_)), CSINTO(mdac_csda_pow_d_into(t_, a, 0.5)),
           da::pow(ra, 0.5));
    S ab;
    OK(mdac_csda_abs(a, ab.out()));
    require_same(ab, da::abs(ra));

    EX sa(parse("a")), sb(parse("b"));
    const mdac_expr* ks[] = {sa, sb};
    const double xs[] = {0.3, -1.25};
    C r;
    OK(mdac_csda_evaluate(a, 2, ks, xs, r.out()));
    SymEngine::vec_basic syms = {SymEngine::symbol("a"), SymEngine::symbol("b")};
    require_same(r, da::evaluate(ra, syms, {0.3, -1.25}), 0.0);
    C r2;
    REQUIRE(mdac_csda_evaluate(a, 1, ks + 1, xs, r2.out()) == MDAC_ERR_VALUE);   // a has no value
    REQUIRE(r2.p == nullptr);
}

TEST_CASE("capi: CSDA lists and cd_composition equal C++", "[capi][sym][csda][algo]") {
    Envs envs(400);
    mdac_env* e = envs.env;
    const RCSDA ra = csref(ca, cb), rb = csref(cc, cd);
    CS a(csfrom(e, ca, cb)), b(csfrom(e, cc, cd));

    SECTION("lists") {
        CSL l;
        OK(mdac_csdalist_new(l.out()));
        REQUIRE(mdac_csdalist_length(l) == 0);
        REQUIRE(mdac_csdalist_env(l) == nullptr);
        OK(mdac_csdalist_push(l, a));
        OK(mdac_csdalist_push(l, b));
        REQUIRE(mdac_csdalist_env(l) == e);
        require_same(l, std::vector<RCSDA>{ra, rb});
        OK(mdac_csdalist_set(l, 1, a));
        require_same(l, std::vector<RCSDA>{ra, ra});
        CS g;
        REQUIRE(mdac_csdalist_get(l, 2, g.out()) == MDAC_ERR_INDEX);
        REQUIRE(mdac_csdalist_set(l, 2, a) == MDAC_ERR_INDEX);
        const mdac_csda* vs[] = {b, a};
        CSL f;
        OK(mdac_csdalist_from(vs, 2, f.out()));
        require_same(f, std::vector<RCSDA>{rb, ra});
    }
    SECTION("cd_composition, three forms") {
        // Arguments: a constant plus a small linear part per variable.
        const std::vector<const char*> v0r = {"7/10", "1/10"}, v0i = {"0", "1/20"};
        const std::vector<const char*> v1r = {"3/10", "0", "1/20"}, v1i = {"1/50"};
        const std::vector<const char*> v2r = {"1/5", "0", "0", "0", "1/10"}, v2i = {"0"};
        const RCSDA rv0 = csref(v0r, v0i), rv1 = csref(v1r, v1i), rv2 = csref(v2r, v2i);
        CS v0(csfrom(e, v0r, v0i)), v1(csfrom(e, v1r, v1i)), v2(csfrom(e, v2r, v2i));
        const mdac_csda* cv[] = {v0, v1, v2};
        const mdac_csda* cm[] = {a, b};
        CSL lv, lm, o1, o2;
        OK(mdac_csdalist_from(cv, 3, lv.out()));
        OK(mdac_csdalist_from(cm, 2, lm.out()));
        std::vector<RCSDA> rv = {rv0, rv1, rv2}, rm = {ra, rb};

        OK(mdac_csdalist_compose(lm, lv, o1.out()));
        std::vector<RCSDA> want(2);
        da::cd_composition(rm, rv, want);
        require_same(o1, want);

        S s0(sfrom(e, ca)), s1(sfrom(e, cc));
        const mdac_sda* sm[] = {s0, s1};
        SL lsm;
        OK(mdac_sdalist_from(sm, 2, lsm.out()));
        std::vector<da::SDA> rsm = {sref(ca), sref(cc)};
        OK(mdac_sdalist_compose_c(lsm, lv, o2.out()));
        std::vector<RCSDA> want2(2);
        da::cd_composition(rsm, rv, want2);
        require_same(o2, want2);

        S w0(sfrom(e, v0r)), w1(sfrom(e, v1r)), w2(sfrom(e, v2r));
        const mdac_sda* sv[] = {w0, w1, w2};
        SL lsv;
        CSL o3;
        OK(mdac_sdalist_from(sv, 3, lsv.out()));
        std::vector<da::SDA> rsv = {sref(v0r), sref(v1r), sref(v2r)};
        OK(mdac_csdalist_compose_s(lm, lsv, o3.out()));
        std::vector<RCSDA> want3(2);
        da::cd_composition(rm, rsv, want3);
        require_same(o3, want3);

        CSL bad, empty, oe;
        OK(mdac_csdalist_from(cv, 2, bad.out()));
        CSL ob;
        REQUIRE(mdac_csdalist_compose(lm, bad, ob.out()) == MDAC_ERR_VALUE);
        REQUIRE(ob.p == nullptr);
        OK(mdac_csdalist_new(empty.out()));
        OK(mdac_csdalist_compose(empty, lv, oe.out()));
        REQUIRE(mdac_csdalist_length(oe) == 0);
    }
}

TEST_CASE("capi: CSDA envs: mixing fails, freeing after mdac_clear is safe", "[capi][sym][csda]") {
    mdac_csda* c = nullptr;
    mdac_csdalist* l = nullptr;
    {
        Envs envs;
        OK(mdac_csda_new_z(envs.env, 1.0, 2.0, &c));
        OK(mdac_csdalist_new(&l));
        OK(mdac_csdalist_push(l, c));
        mdac_env* other = nullptr;
        OK(mdac_env_make(kOrder, kNvars, 50, 0, &other));
        CS t;
        OK(mdac_csda_new_z(other, 1.0, 0.0, t.out()));
        S s;
        OK(mdac_sda_new_d(other, 2.0, s.out()));
        S s2;
        OK(mdac_sda_new_d(envs.env, 2.0, s2.out()));
        CS r;
        REQUIRE(mdac_csda_add(c, t, r.out()) == MDAC_ERR_ENV);
        REQUIRE(mdac_csda_mul_s(c, s, r.out()) == MDAC_ERR_ENV);
        REQUIRE(mdac_csda_add_into(t, c, c) == MDAC_ERR_ENV);
        REQUIRE(mdac_csda_new(s2, s, r.out()) == MDAC_ERR_ENV);
        REQUIRE(mdac_csda_set_real(c, s) == MDAC_ERR_ENV);
        REQUIRE(mdac_csdalist_push(l, t) == MDAC_ERR_ENV);
        mdac_csda_free(t.p);
        t.p = nullptr;
        mdac_sda_free(s.p);
        s.p = nullptr;
        OK(mdac_env_close(other));
    }
    CS r;
    REQUIRE(mdac_csda_exp(c, r.out()) == MDAC_ERR_ENV);
    mdac_csda_free(c);
    mdac_csdalist_free(l);
}

TEST_CASE("capi: import of SDA and CSDA, promote_to", "[capi][sym][multienv]") {
    Envs envs;
    mdac_env *b = nullptr, *small = nullptr;
    OK(mdac_env_make(kOrder, kNvars, 20, 0, &b));
    OK(mdac_env_make(2, 2, 20, 0, &small));
    H x(from(envs.env, input(1.5, 0.7)));
    S p, ref;
    OK(mdac_nda_promote_to(b, x, p.out()));
    OK(mdac_sda_promote(x, ref.out()));
    REQUIRE(mdac_sda_env(p) == b);
    mdac_env* cur = nullptr;
    OK(mdac_env_current(&cur));
    REQUIRE(cur == envs.env);
    REQUIRE(scoeffs(p) == scoeffs(ref));
    S q;
    OK(mdac_sda_import(envs.env, p, q.out()));
    REQUIRE(mdac_sda_env(q) == envs.env);
    REQUIRE(scoeffs(q) == scoeffs(ref));

    CS c, cb;
    OK(mdac_csda_new(q, ref, c.out()));
    OK(mdac_csda_import(b, c, cb.out()));
    REQUIRE(mdac_csda_env(cb) == b);
    S re, im;
    OK(mdac_csda_real(cb, re.out()));
    OK(mdac_csda_imag(cb, im.out()));
    REQUIRE(scoeffs(re) == scoeffs(ref));
    REQUIRE(scoeffs(im) == scoeffs(ref));

    S bad;
    REQUIRE(mdac_nda_promote_to(small, x, bad.out()) == MDAC_ERR_VALUE);
    REQUIRE(mdac_sda_import(small, q, bad.out()) == MDAC_ERR_VALUE);
    CS cbad;
    REQUIRE(mdac_csda_import(small, c, cbad.out()) == MDAC_ERR_VALUE);
    mdac_sda_free(p.p);
    p.p = nullptr;
    mdac_sda_free(re.p);
    re.p = nullptr;
    mdac_sda_free(im.p);
    im.p = nullptr;
    mdac_csda_free(cb.p);
    cb.p = nullptr;
    OK(mdac_env_close(b));
    REQUIRE(mdac_nda_promote_to(b, x, bad.out()) == MDAC_ERR_ENV);
    REQUIRE(mdac_sda_import(b, q, bad.out()) == MDAC_ERR_ENV);
    REQUIRE(mdac_csda_import(b, c, cbad.out()) == MDAC_ERR_ENV);
    OK(mdac_env_close(small));
}

#endif // DA_WITH_SYMBOLIC
