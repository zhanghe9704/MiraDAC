/**
 * @file example_complex_symbolic.cc
 * @brief Example: symbolic complex DA (std::complex<SDA> = CSDA).
 *
 * @details Demonstrates:
 *   - Constructing a CSDA map with symbolic coefficients.
 *   - Complex DA arithmetic and transcendental functions.
 *   - Separating real/imaginary parts.
 *   - Evaluating the symbolic result to a CNDA at concrete numeric values.
 *
 * CSDA is the supported representation for complex symbolic DA: the
 * imaginary unit lives in the std::complex wrapper, so real and imaginary
 * parts stay separable and evaluate() can return a CNDA. Do not put
 * SymEngine::I inside an SDA coefficient — see the CSDA alias in da/da.h.
 *
 * Build (symbolic): cmake ... -DWITH_SYMBOLIC=ON && cmake --build ...
 *
 * Stage C6 per DEVELOPMENT_PLAN_COMPLEX_SYMBOLIC.md.
 */

#ifdef DA_WITH_SYMBOLIC

#include "da/da.h"

#include <symengine/expression.h>
#include <symengine/symbol.h>

#include <complex>
#include <iostream>

using namespace da;
namespace SE = SymEngine;

int main() {
    // Initialize: order 3, 2 variables, pool 300
    da::da_init(3, 2, 300);

    // Symbolic parameters for the linear DA coefficients
    SE::Expression a = SE::Expression(SE::symbol("a"));
    SE::Expression b = SE::Expression(SE::symbol("b"));

    // Build a symbolic complex DA map:
    //   z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2)
    // The constant parts (1.5, 0.5) are numeric so that scalar DA functions
    // (exp, sqrt, etc.) can branch on them; higher-order coefficients are symbolic.
    SDA re = SDA(SE::Expression(1.5))
             + a * promote(da::base[0])
             + SDA(SE::Expression(0.1)) * promote(da::base[1]);
    SDA im = SDA(SE::Expression(0.5))
             + b * promote(da::base[0])
             + SDA(SE::Expression(-0.05)) * promote(da::base[1]);
    CSDA z(re, im);

    std::cout << "=== Symbolic Complex DA example ===\n\n";
    std::cout << "z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2)\n\n";

    // Arithmetic
    CSDA z2 = z * z;
    std::cout << "z^2 real constant part: "
              << static_cast<double>(get_real(z2).con()) << "\n";  // 1.5^2 - 0.5^2 = 2.0
    std::cout << "z^2 imag constant part: "
              << static_cast<double>(get_imag(z2).con()) << "\n";  // 2*1.5*0.5 = 1.5

    // Transcendental function
    CSDA ez = da::exp(z);
    std::cout << "\nexp(z) real constant part: "
              << static_cast<double>(get_real(ez).con())  // exp(1.5)*cos(0.5)
              << " (expected ~" << std::exp(1.5) * std::cos(0.5) << ")\n";
    std::cout << "exp(z) imag constant part: "
              << static_cast<double>(get_imag(ez).con())  // exp(1.5)*sin(0.5)
              << " (expected ~" << std::exp(1.5) * std::sin(0.5) << ")\n";

    // Separate real and imaginary parts
    const SDA& re_ez = get_real(ez);
    const SDA& im_ez = get_imag(ez);
    std::cout << "\nReal part of exp(z) at order 0: " << static_cast<double>(re_ez.con()) << "\n";
    std::cout << "Imag part of exp(z) at order 0: " << static_cast<double>(im_ez.con()) << "\n";

    // Evaluate at concrete values: a = 0.3, b = 0.2
    SE::map_basic_basic vals;
    vals[SE::symbol("a")] = SE::real_double(0.3);
    vals[SE::symbol("b")] = SE::real_double(0.2);

    CNDA ez_num = evaluate(ez, vals);

    // Compare with direct numerical computation
    NDA re_n = NDA(1.5) + NDA(0.3) * da::base[0] + NDA(0.1) * da::base[1];
    NDA im_n = NDA(0.5) + NDA(0.2) * da::base[0] + NDA(-0.05) * da::base[1];
    CNDA n(re_n, im_n);
    CNDA ez_direct = da::exp(n);

    bool ok = da::compare_cd_vectors(ez_num, ez_direct, 1e-12);
    std::cout << "\nEvaluation cross-check (exp): " << (ok ? "PASS" : "FAIL") << "\n";

    // sqrt example
    CSDA sz = da::sqrt(z);
    CNDA sz_num = evaluate(sz, vals);
    CNDA sz_direct = da::sqrt(n);
    bool ok2 = da::compare_cd_vectors(sz_num, sz_direct, 1e-12);
    std::cout << "Evaluation cross-check (sqrt): " << (ok2 ? "PASS" : "FAIL") << "\n";

    std::cout << "\nAll CSDA operations completed successfully.\n";
    return 0;
}

#else

#include <iostream>
int main() {
    std::cout << "This example requires WITH_SYMBOLIC=ON\n";
    return 0;
}

#endif // DA_WITH_SYMBOLIC
