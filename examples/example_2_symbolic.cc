/**
 * @file example_2_symbolic.cc
 * @brief Symbolic DA: plug in values using eval().
 *
 * Ported from ref/tpsa_sym/examples/example_2_eval.cc.
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
    int da_order = 3;
    int da_dim   = 3;
    int da_pool  = 400;
    da::da_init(static_cast<unsigned>(da_order),
                static_cast<unsigned>(da_dim),
                static_cast<unsigned>(da_pool),
                /*table=*/true);

    Expression x("x0"), y("y0"), z("z0");
    auto& b = da::base;

    SDA da1 = Expression(1) + (Expression(1)+x)*b[0] + y*b[1] + (z-Expression(0.5))*b[2];
    cout << da1;

    // Evaluate: x=0.1, y=0.2, z=0.3
    SymEngine::map_basic_basic subs_map;
    subs_map[x] = Expression("0.1");
    subs_map[y] = Expression("0.2");
    subs_map[z] = Expression("0.3");
    da1.eval(subs_map);

    cout << da1;
    return 0;
}

#else
#include <iostream>
int main() {
    std::cout << "example_2_symbolic: built without DA_WITH_SYMBOLIC, skipping." << std::endl;
    return 0;
}
#endif
