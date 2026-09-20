/**
 * @file example_1_symbolic.cc
 * @brief Symbolic DA fundamental operations.
 *
 * Ported from ref/tpsa_sym/examples/example_1_fundamental.cc.
 * Requires DA_WITH_SYMBOLIC.
 */

#ifdef DA_WITH_SYMBOLIC

#include "da/da.h"
#include <iostream>

using SymEngine::Expression;
using da::SDA;
using std::cout;
using std::endl;

int main() {
    // Initialize DA environment
    int da_order = 4;
    int da_dim   = 3;
    int da_pool  = 400;
    da::da_init(static_cast<unsigned>(da_order),
                static_cast<unsigned>(da_dim),
                static_cast<unsigned>(da_pool),
                /*table=*/true);

    // Symbolic variables
    Expression x("x"), y("y"), z("z");

    // Bases: da::base[0], da::base[1], da::base[2]
    auto& b = da::base;

    // Construct symbolic DA vectors
    SDA da1 = Expression(1) + (Expression(1)+x)*b[0] + y*b[1] + (z-Expression(0.5))*b[2];
    cout << "Symbolic DA vector 1:" << endl << da1;

    SDA da2 = Expression(3.3) + (Expression(0.5)+x)*b[0] + y*y*b[1]
            + (x+z+Expression(1.1))*b[2];
    cout << "Symbolic DA vector 2:" << endl << da2;

    cout << "Summation: vec 1 + vec 2:" << endl;
    SDA da3 = da1 + da2;
    cout << da3;

    cout << "Multiplication: vec 1 * vec 2:" << endl;
    da3 = da1 * da2;
    cout << da3;

    cout << "Inverse of vec 1:" << endl;
    da3 = Expression(1) / da1;
    cout << da3;

    cout << "Square root of vec 1:" << endl;
    da3 = da::sqrt(da1);
    cout << da3;

    return 0;
}

#else // DA_WITH_SYMBOLIC not defined
#include <iostream>
int main() {
    std::cout << "example_1_symbolic: built without DA_WITH_SYMBOLIC, skipping." << std::endl;
    return 0;
}
#endif
