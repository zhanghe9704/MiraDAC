// test_capi.cc — the C API against direct C++ results (plan T1.2-T1.5).
//
// The C API library holds its own hidden copy of the da library, so the C++
// reference vectors below live in a separate default env of the same shape.
#include "catch.hpp"
#include "miradac.h"

#include "da/da.h"

#include <cmath>
#include <complex>
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
