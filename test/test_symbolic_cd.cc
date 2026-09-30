/**
 * @file test_symbolic_cd.cc
 * @brief Tests for symbolic complex DA (std::complex<SDA> = CSDA).
 *
 * @details Stage C5: correctness validated by evaluation cross-check.
 *   For each operation:
 *     1. Build symbolic CSDA inputs with NUMERIC constant parts and symbolic
 *        higher-order coefficients (e.g., c0 = 0.5, c1 = symbol "b" * x1).
 *        This ensures the scalar SDA functions can branch on numeric c0 while
 *        the symbolic higher-order coefficients exercise the symbolic path.
 *     2. Compute the symbolic result Rs = op(inputs).
 *     3. Substitute concrete numeric values for the symbols.
 *     4. Evaluate to CNDA R_from_sym.
 *     5. Compute the same op numerically on CNDA inputs R_num.
 *     6. REQUIRE coefficient-wise agreement within 1e-10.
 *
 * The entire file is guarded by DA_WITH_SYMBOLIC so the OFF build is unaffected.
 */

#ifdef DA_WITH_SYMBOLIC

#include "catch.hpp"
#include <sstream>
#include "da/da.h"

#include <symengine/expression.h>
#include <symengine/symbol.h>

#include <complex>
#include <cmath>
#include <vector>
#include <string>

using namespace da;
namespace SE = SymEngine;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static SE::Expression make_sym(const std::string& name) {
    return SE::Expression(SE::symbol(name));
}

// Build a symbol map: symbol_name -> double
static SE::map_basic_basic make_vals(
    std::initializer_list<std::pair<std::string, double>> vals) {
    SE::map_basic_basic m;
    for (auto& kv : vals)
        m[SE::symbol(kv.first)] = SE::real_double(kv.second);
    return m;
}

// Build a CSDA with numeric constant part and symbolic linear coefficients:
//   re = c0r + b1r*x1 + b2r*x2,  im = c0i + b1i*x1 + b2i*x2
// c0r, c0i are numeric; b1r, b2r, b1i, b2i are SymEngine Expression (may be symbols).
static CSDA make_csda(
    double c0r, const SE::Expression& b1r, const SE::Expression& b2r,
    double c0i, const SE::Expression& b1i, const SE::Expression& b2i)
{
    SDA re = SDA(SE::Expression(c0r))
             + b1r * promote(da::base[0])
             + b2r * promote(da::base[1]);
    SDA im = SDA(SE::Expression(c0i))
             + b1i * promote(da::base[0])
             + b2i * promote(da::base[1]);
    return CSDA(re, im);
}

// Build corresponding CNDA at given numeric values for the symbols
static CNDA make_cnda(
    double c0r, double b1r_v, double b2r_v,
    double c0i, double b1i_v, double b2i_v)
{
    NDA re = NDA(c0r) + NDA(b1r_v) * da::base[0] + NDA(b2r_v) * da::base[1];
    NDA im = NDA(c0i) + NDA(b1i_v) * da::base[0] + NDA(b2i_v) * da::base[1];
    return CNDA(re, im);
}

// Compare two CNDA coefficient-by-coefficient within eps
static bool cnda_eq(CNDA& a, CNDA& b, double eps = 1e-10) {
    return compare_cd_vectors(a, b, eps);
}

// Evaluate CSDA at given symbol values -> CNDA
static CNDA eval_csda(const CSDA& z, const SE::map_basic_basic& vals) {
    return evaluate(z, vals);
}

// ---------------------------------------------------------------------------
// TEST CASES
// ---------------------------------------------------------------------------

TEST_CASE("Symbolic CD: structural check and aliases", "[symbolic_cd]") {
    da::da_init(3, 2, 200);

    // CNDA <-> CSDA type aliases exist
    static_assert(std::is_same_v<CNDA, std::complex<NDA>>,  "CNDA alias wrong");
    static_assert(std::is_same_v<CSDA, std::complex<SDA>>,  "CSDA alias wrong");

    SE::Expression b = make_sym("b");
    CSDA z = make_csda(1.0, b, SE::Expression(0.0), 0.5, SE::Expression(0.2), SE::Expression(0.0));

    // get_real/get_imag roundtrip
    const SDA& re = get_real(z);
    const SDA& im = get_imag(z);

    // Constant parts are numeric
    REQUIRE(static_cast<double>(re.con()) == Approx(1.0));
    REQUIRE(static_cast<double>(im.con()) == Approx(0.5));

    // Arithmetic compiles and produces correctly-zero diff
    CSDA diff = z - z;
    REQUIRE(static_cast<double>(diff.real().con()) == Approx(0.0));
    REQUIRE(static_cast<double>(diff.imag().con()) == Approx(0.0));
}

TEST_CASE("Symbolic CD: evaluation cross-check — arithmetic operators", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    SE::Expression b1 = make_sym("b1"), b2 = make_sym("b2");
    SE::Expression c1 = make_sym("c1"), c2 = make_sym("c2");
    SE::Expression e1 = make_sym("e1"), e2 = make_sym("e2");
    SE::Expression f1 = make_sym("f1"), f2 = make_sym("f2");

    // z1 = (0.8 + b1*x1 + b2*x2) + i*(0.3 + c1*x1 + c2*x2)
    CSDA z1 = make_csda(0.8, b1, b2, 0.3, c1, c2);
    // z2 = (0.5 + e1*x1 + e2*x2) + i*(0.6 + f1*x1 + f2*x2)
    CSDA z2 = make_csda(0.5, e1, e2, 0.6, f1, f2);

    // Concrete values for the symbols
    double b1_v = 0.3, b2_v = -0.1, c1_v = 0.2, c2_v = 0.4;
    double e1_v = -0.15, e2_v = 0.05, f1_v = 0.1, f2_v = -0.2;

    auto vals = make_vals({{"b1", b1_v}, {"b2", b2_v}, {"c1", c1_v}, {"c2", c2_v},
                           {"e1", e1_v}, {"e2", e2_v}, {"f1", f1_v}, {"f2", f2_v}});

    CNDA n1 = make_cnda(0.8, b1_v, b2_v, 0.3, c1_v, c2_v);
    CNDA n2 = make_cnda(0.5, e1_v, e2_v, 0.6, f1_v, f2_v);

    double eps = 1e-10;

    // Addition
    {
        CSDA rs = z1 + z2;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 + n2;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // Subtraction
    {
        CSDA rs = z1 - z2;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 - n2;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // Multiplication
    {
        CSDA rs = z1 * z2;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 * n2;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // Division
    {
        CSDA rs = z1 / z2;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 / n2;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // z + complex<double>
    {
        std::complex<double> c(0.5, -0.3);
        CSDA rs = z1 + c;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 + c;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // z * double scalar
    {
        CSDA rs = z1 * 2.5;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = n1 * 2.5;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // double / z
    {
        CSDA rs = 1.0 / z2;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = 1.0 / n2;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // (-1)*z negation
    {
        CSDA rs = (-1.0) * z1;
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = (-1.0) * n1;
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
}

TEST_CASE("Symbolic CD: evaluation cross-check — exp, sqrt, log", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    SE::Expression b1 = make_sym("b1"), b2 = make_sym("b2");
    SE::Expression c1 = make_sym("c1"), c2 = make_sym("c2");

    // z = (0.5 + b1*x1 + b2*x2) + i*(0.3 + c1*x1 + c2*x2)
    CSDA z = make_csda(0.5, b1, b2, 0.3, c1, c2);

    double b1_v = 0.10, b2_v = 0.05, c1_v = 0.08, c2_v = -0.04;
    auto vals = make_vals({{"b1", b1_v}, {"b2", b2_v}, {"c1", c1_v}, {"c2", c2_v}});
    CNDA n = make_cnda(0.5, b1_v, b2_v, 0.3, c1_v, c2_v);

    double eps = 1e-9;

    // exp
    {
        CSDA rs = da::exp(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::exp(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // sqrt
    {
        CSDA rs = da::sqrt(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::sqrt(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // log
    {
        CSDA rs = da::log(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::log(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
}

TEST_CASE("Symbolic CD: evaluation cross-check — asin, acos, atan", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");

    // z = (0.3 + b1*x1) + i*(0.2 + c1*x1), both small to keep |z| < 1
    CSDA z = make_csda(0.3, b1, SE::Expression(0.0), 0.2, c1, SE::Expression(0.0));
    double b1_v = 0.05, c1_v = 0.03;
    auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
    CNDA n = make_cnda(0.3, b1_v, 0.0, 0.2, c1_v, 0.0);

    double eps = 1e-9;

    // asin
    {
        CSDA rs = da::asin(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::asin(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // acos
    {
        CSDA rs = da::acos(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::acos(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // atan
    {
        CSDA rs = da::atan(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::atan(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
}

TEST_CASE("Symbolic CD: evaluation cross-check — asinh, acosh, atanh", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");

    double eps = 1e-9;

    // asinh, atanh: small constant part
    {
        CSDA z = make_csda(0.5, b1, SE::Expression(0.0), 0.3, c1, SE::Expression(0.0));
        double b1_v = 0.08, c1_v = 0.05;
        auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
        CNDA n = make_cnda(0.5, b1_v, 0.0, 0.3, c1_v, 0.0);

        CSDA rs = da::asinh(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::asinh(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));

        CSDA rs2 = da::atanh(z);
        CNDA r_sym2 = eval_csda(rs2, vals);
        CNDA r_num2 = da::atanh(n);
        REQUIRE(cnda_eq(r_sym2, r_num2, eps));
    }
    // acosh: constant part > 1
    {
        CSDA z = make_csda(1.5, b1, SE::Expression(0.0), 0.2, c1, SE::Expression(0.0));
        double b1_v = 0.1, c1_v = 0.05;
        auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
        CNDA n = make_cnda(1.5, b1_v, 0.0, 0.2, c1_v, 0.0);

        CSDA rs = da::acosh(z);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::acosh(n);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
}

TEST_CASE("Symbolic CD: evaluation cross-check — pow", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");
    CSDA z = make_csda(0.5, b1, SE::Expression(0.0), 0.3, c1, SE::Expression(0.0));
    double b1_v = 0.08, c1_v = 0.05;
    auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
    CNDA n = make_cnda(0.5, b1_v, 0.0, 0.3, c1_v, 0.0);

    double eps = 1e-9;

    // pow(z, int)
    {
        CSDA rs = da::pow(z, 3);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::pow(n, 3);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
    // pow(z, double)
    {
        CSDA rs = da::pow(z, 2.5);
        CNDA r_sym = eval_csda(rs, vals);
        CNDA r_num = da::pow(n, 2.5);
        REQUIRE(cnda_eq(r_sym, r_num, eps));
    }
}

TEST_CASE("Symbolic CD: abs (symbolic magnitude)", "[symbolic_cd]") {
    da::da_init(3, 2, 300);

    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");
    CSDA z = make_csda(0.6, b1, SE::Expression(0.0), 0.4, c1, SE::Expression(0.0));
    double b1_v = 0.1, c1_v = 0.08;
    auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
    CNDA n = make_cnda(0.6, b1_v, 0.0, 0.4, c1_v, 0.0);

    // Symbolic abs returns SDA: sqrt(re^2 + im^2)
    SDA rs = da::abs(z);
    NDA r_sym_nda = evaluate(rs, vals);

    // Numerical abs: also sqrt(re^2+im^2) — compute numerically
    NDA r_num_mag = da::sqrt(get_real(n) * get_real(n) + get_imag(n) * get_imag(n));

    REQUIRE(compare_da_vectors(r_sym_nda, r_num_mag, 1e-10));
}

TEST_CASE("Symbolic CD: promote/evaluate interop", "[symbolic_cd]") {
    da::da_init(3, 2, 300);

    // promote CNDA -> CSDA -> evaluate back
    CNDA orig(NDA(1.5) + NDA(0.3) * da::base[0] + NDA(-0.1) * da::base[1],
              NDA(0.7) + NDA(0.2) * da::base[0] + NDA(-0.05) * da::base[1]);

    CSDA promoted = promote(orig);

    SE::map_basic_basic empty_map;
    CNDA roundtrip = evaluate(promoted, empty_map);
    REQUIRE(cnda_eq(orig, roundtrip, 1e-14));
}

TEST_CASE("Symbolic CD: cd_composition cross-check", "[symbolic_cd]") {
    da::da_init(3, 2, 800);

    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");

    // Symbolic complex map: z[0] = (1.2 + b1*x1) + i*(0.4 + c1*x1)
    //                        z[1] = 0.5 + i*0.2 (constant)
    CSDA z0 = make_csda(1.2, b1, SE::Expression(0.0), 0.4, c1, SE::Expression(0.0));
    CSDA z1 = make_csda(0.5, SE::Expression(0.0), SE::Expression(0.0),
                         0.2, SE::Expression(0.0), SE::Expression(0.0));
    std::vector<CSDA> sym_in = { z0, z1 };

    // Substitution map: numeric SDA -> CSDA
    NDA v0_re = NDA(0.7) + NDA(0.1) * da::base[0];
    NDA v0_im = NDA(0.05) * da::base[0];
    NDA v1_re = NDA(0.3) + NDA(0.05) * da::base[1];
    NDA v1_im = NDA(0.02) * da::base[1];
    std::vector<CSDA> cvec = {
        CSDA(promote(v0_re), promote(v0_im)),
        CSDA(promote(v1_re), promote(v1_im))
    };

    std::vector<CSDA> sym_out(2);
    da::cd_composition(sym_in, cvec, sym_out);

    // Evaluate symbolically
    double b1_v = 0.3, c1_v = 0.15;
    auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
    CNDA sym_out0_eval = evaluate(sym_out[0], vals);
    CNDA sym_out1_eval = evaluate(sym_out[1], vals);

    // Numerical version
    CNDA n0(NDA(1.2) + NDA(b1_v) * da::base[0], NDA(0.4) + NDA(c1_v) * da::base[0]);
    CNDA n1(NDA(0.5), NDA(0.2));
    std::vector<CNDA> num_in = { n0, n1 };
    std::vector<CNDA> ncvec = {
        CNDA(v0_re, v0_im), CNDA(v1_re, v1_im)
    };
    std::vector<CNDA> num_out(2);
    da::cd_composition(num_in, ncvec, num_out);

    double eps = 1e-9;
    REQUIRE(cnda_eq(sym_out0_eval, num_out[0], eps));
    REQUIRE(cnda_eq(sym_out1_eval, num_out[1], eps));
}

TEST_CASE("Symbolic CD: cd_composition vector<SDA> into complex map", "[symbolic_cd]") {
    da::da_init(3, 2, 500);

    // Test the vector<SDA> into vector<CSDA> overload
    SE::Expression b1 = make_sym("b1"), c1 = make_sym("c1");

    SDA re_sym = SDA(SE::Expression(1.1)) + b1 * promote(da::base[0]);
    SDA im_sym = SDA(SE::Expression(0.5)) + c1 * promote(da::base[1]);
    std::vector<SDA> sym_ivec = { re_sym, im_sym };

    // Numeric composition map
    NDA cv0re = NDA(0.5) + NDA(0.1) * da::base[0];
    NDA cv1re = NDA(0.3);
    NDA cv0im = NDA(0.2);
    NDA cv1im = NDA(0.05) * da::base[1];

    std::vector<CSDA> cmap = {
        CSDA(promote(cv0re), promote(cv0im)),
        CSDA(promote(cv1re), promote(cv1im))
    };
    std::vector<CSDA> sout(2);
    da::cd_composition(sym_ivec, cmap, sout);

    // Also numeric version
    double b1_v = 0.3, c1_v = 0.2;
    NDA re_num = NDA(1.1) + NDA(b1_v) * da::base[0];
    NDA im_num = NDA(0.5) + NDA(c1_v) * da::base[1];
    std::vector<NDA> num_ivec = { re_num, im_num };
    std::vector<CNDA> cmap_num = {
        CNDA(cv0re, cv0im), CNDA(cv1re, cv1im)
    };
    std::vector<CNDA> num_out(2);
    da::cd_composition(num_ivec, cmap_num, num_out);

    auto vals = make_vals({{"b1", b1_v}, {"c1", c1_v}});
    CNDA sout0_eval = evaluate(sout[0], vals);
    CNDA sout1_eval = evaluate(sout[1], vals);

    double eps = 1e-9;
    REQUIRE(cnda_eq(sout0_eval, num_out[0], eps));
    REQUIRE(cnda_eq(sout1_eval, num_out[1], eps));
}

// ---------------------------------------------------------------------------
// CSDA mixed with SDA and with a symbolic scalar (Expression)
// ---------------------------------------------------------------------------
TEST_CASE("Symbolic CD: CSDA with SDA and Expression operands", "[symbolic_cd]") {
    da::da_init(3, 2, 800);

    SE::Expression b1 = make_sym("b1"), b2 = make_sym("b2");
    SE::Expression c1 = make_sym("c1"), c2 = make_sym("c2");
    SE::Expression e1 = make_sym("e1"), k = make_sym("k");
    CSDA z = make_csda(0.8, b1, b2, 0.3, c1, c2);
    SDA s = SDA(SE::Expression(0.5)) + e1 * promote(da::base[0]);

    auto vals = make_vals({{"b1", 0.3}, {"b2", -0.1}, {"c1", 0.2}, {"c2", 0.4},
                           {"e1", -0.15}, {"k", 1.7}});
    CNDA zn = make_cnda(0.8, 0.3, -0.1, 0.3, 0.2, 0.4);
    NDA sn = NDA(0.5) + NDA(-0.15) * da::base[0];
    const double kn = 1.7;

    auto check = [&](const CSDA& sym, CNDA num) {
        CNDA got = eval_csda(sym, vals);
        REQUIRE(cnda_eq(got, num, 1e-12));
    };

    SECTION("CSDA op SDA") {
        check(z + s, zn + sn);  check(s + z, sn + zn);
        check(z - s, zn - sn);  check(s - z, sn - zn);
        check(z * s, zn * sn);  check(s * z, sn * zn);
        check(z / s, zn / sn);  check(s / z, sn / zn);
    }
    SECTION("CSDA op Expression") {
        check(z + k, zn + kn);  check(k + z, kn + zn);
        check(z - k, zn - kn);  check(k - z, kn - zn);
        check(z * k, zn * kn);  check(k * z, kn * zn);
        check(z / k, zn / kn);  check(k / z, kn / zn);
    }
    SECTION("operator<< prints both parts") {
        std::ostringstream os;
        os << z;
        const std::string out = os.str();
        REQUIRE(out.find("Real part") != std::string::npos);
        REQUIRE(out.find("Imaginary part") != std::string::npos);
        REQUIRE(out.find("b1") != std::string::npos);
        REQUIRE(out.find("c1") != std::string::npos);
    }
}

#endif // DA_WITH_SYMBOLIC
