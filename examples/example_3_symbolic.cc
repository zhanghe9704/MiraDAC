/**
 * @file example_3_symbolic.cc
 * @brief Symbolic DA: eval_funs for efficient repeated evaluation.
 *
 * Ported from ref/tpsa_sym/examples/example_3_eval_funs.cc.
 * Requires DA_WITH_SYMBOLIC.
 */

#ifdef DA_WITH_SYMBOLIC

#include "da/da.h"
#include <array>
#include <iostream>
#include <vector>

using SymEngine::Expression;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::symbol;
using da::SDA;
using std::cout;
using std::endl;
using std::vector;
using std::array;

int main() {
    int da_order = 3;
    int da_dim   = 3;
    int da_pool  = 400;
    da::da_init(static_cast<unsigned>(da_order),
                static_cast<unsigned>(da_dim),
                static_cast<unsigned>(da_pool),
                /*table=*/true);

    RCP<const Basic> x = symbol("x0");
    RCP<const Basic> y = symbol("y0");
    RCP<const Basic> z = symbol("z0");

    auto& b = da::base;

    SDA da1 = Expression(1) + (Expression(1)+x)*b[0] + Expression(y)*b[1]
            + (Expression(z)-Expression(0.5))*b[2];
    cout << da1;

    // Per-coefficient visitors
    vector<RCP<const Basic>> vars{x, y, z};
    vector<SymEngine::LambdaRealDoubleVisitor> vec;
    da1.eval_funs(vars, vec);
    cout << "Per-coefficient evaluation at x=0.1, y=0.2, z=0.3:" << endl;
    for (auto& v : vec) {
        cout << v.call({0.1, 0.2, 0.3}) << endl;
    }

    // Single vectorised visitor
    SymEngine::LambdaRealDoubleVisitor v;
    da1.eval_funs(vars, v);
    vector<double> results(static_cast<std::size_t>(da1.full_length()));
    array<double,3> inputs{0.1, 0.2, 0.3};
    cout << "==============" << endl;
    v.call(&results[0], &inputs[0]);
    for (auto r : results)
        cout << r << endl;

    return 0;
}

#else
#include <iostream>
int main() {
    std::cout << "example_3_symbolic: built without DA_WITH_SYMBOLIC, skipping." << std::endl;
    return 0;
}
#endif
