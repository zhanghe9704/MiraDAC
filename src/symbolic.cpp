/**
 * @file symbolic.cpp
 * @brief Symbolic utility functions (is_zero, simplified_expr).
 *
 * @details Ported from ref/tpsa_sym/src/symbolic.cc.
 * Compiled only when DA_WITH_SYMBOLIC is defined.
 */

#ifdef DA_WITH_SYMBOLIC

#include "da/symbolic_ops.h"

#include <symengine/expression.h>
#include <symengine/simplify.h>

using SymEngine::Expression;
using SymEngine::simplify;
using SymEngine::expand;

namespace da {

bool is_zero(Expression x) {
    Expression y = simplify(expand(x));
    return y == Expression(0) || y == Expression(0.0);
}

void simplified_expr(Expression& x) {
    x = simplify(expand(x));
}

} // namespace da

#endif // DA_WITH_SYMBOLIC
