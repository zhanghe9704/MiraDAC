/**
 * @file example_interop.cc
 * @brief Interoperability example: NDA * SDA -> SDA, then evaluate -> NDA.
 *
 * Demonstrates Stage 7 interop:
 *   1. Build an NDA (numerical DA vector).
 *   2. Build an SDA (symbolic DA vector with a symbolic coefficient).
 *   3. Multiply them: NDA * SDA -> SDA (via automatic promote()).
 *   4. Evaluate the result back to NDA by substituting the symbol value.
 *
 * Build: cmake -DWITH_SYMBOLIC=ON ... && cmake --build ... --target example_interop
 * Run:   LD_LIBRARY_PATH=<symengine_lib> ./examples/example_interop
 */

#ifdef DA_WITH_SYMBOLIC

#include "da/da.h"
#include <iostream>
#include <iomanip>

int main() {
    using namespace da;
    using SymEngine::Expression;

    // Initialise: order=4, 2 variables, pool=500, table=true (needed for symbolic)
    da_init(4, 2, 500, true);

    std::cout << "=== MiraDAC Interop Example ===" << std::endl;
    std::cout << "order=" << NDA::order() << "  nv=" << NDA::dim()
              << "  full_len=" << NDA::full_length() << std::endl << std::endl;

    // --- 1. Numerical DA vector ---
    // nx = x (first basis variable)
    NDA nx = base[0];
    std::cout << "NDA: nx = base[0]" << std::endl;
    std::cout << nx << std::endl;

    // --- 2. Symbolic DA vector ---
    // Create a symbolic coefficient 'alpha'
    Expression alpha("alpha");

    // sy = alpha * base[1]  (alpha times the second basis variable)
    SDA sy = alpha * base[1];
    std::cout << "\nSDA: sy = alpha * base[1]" << std::endl;
    std::cout << sy << std::endl;

    // --- 3. Mixed multiplication: NDA * SDA -> SDA ---
    // result_sda = nx * sy = x * (alpha*y) = alpha*x*y
    SDA result_sda = nx * sy;
    std::cout << "\nSDA result = nx * sy  (NDA * SDA -> SDA, via promote):" << std::endl;
    std::cout << result_sda << std::endl;

    // --- 4. Evaluate: substitute alpha=3.0, get NDA ---
    SymEngine::map_basic_basic values = {
        { alpha.get_basic(), Expression(3.0).get_basic() }
    };
    NDA result_nda = evaluate(result_sda, values);
    std::cout << "\nNDA result = evaluate(result_sda, {alpha=3.0}):" << std::endl;
    std::cout << result_nda << std::endl;

    // Verify: coefficient of x*y (index {1,1}) should be 3.0
    std::vector<int> idx{1, 1};
    double coeff = result_nda.element(idx);
    std::cout << "\nCoefficient of x*y = " << std::fixed << std::setprecision(6)
              << coeff << "  (expected 3.0)" << std::endl;

    if (std::abs(coeff - 3.0) < 1e-12) {
        std::cout << "PASSED: coefficient matches 3.0" << std::endl;
    } else {
        std::cout << "FAILED: coefficient mismatch!" << std::endl;
        da_clear();
        return 1;
    }

    da_clear();
    return 0;
}

#else

#include <iostream>
int main() {
    std::cout << "example_interop requires DA_WITH_SYMBOLIC=ON" << std::endl;
    return 0;
}

#endif // DA_WITH_SYMBOLIC
