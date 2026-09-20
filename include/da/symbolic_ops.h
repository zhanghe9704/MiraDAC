/**
 * @file symbolic_ops.h
 * @brief Scalar operator overloads for SymEngine::Expression.
 *
 * @details These overloads allow int, unsigned int, and double to interoperate
 *   with SymEngine::Expression using +, -, *, /.  They are needed so that the
 *   generic engine kernels (which compute e.g. T{1}/c) compile correctly when
 *   T = SymEngine::Expression.
 *
 * Ported from ref/tpsa_sym/include/symbolic.h.
 * Guarded by DA_WITH_SYMBOLIC.
 */

#pragma once

#ifdef DA_WITH_SYMBOLIC

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/simplify.h>

namespace da {

// -----------------------------------------------------------------------
// is_zero: symbolic zero check using simplify + expand
// -----------------------------------------------------------------------
bool is_zero(SymEngine::Expression x);

// -----------------------------------------------------------------------
// simplified_expr: simplify in-place
// -----------------------------------------------------------------------
void simplified_expr(SymEngine::Expression& x);

} // namespace da

// -----------------------------------------------------------------------
// Scalar operator overloads for SymEngine::Expression
// (in global namespace, matching the reference pattern)
// -----------------------------------------------------------------------

// --- operator/ ---
inline SymEngine::Expression operator/(int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) / expr;
}
inline SymEngine::Expression operator/(const SymEngine::Expression expr, int i) {
    return expr / SymEngine::Expression(i);
}
inline SymEngine::Expression operator/(unsigned int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) / expr;
}
inline SymEngine::Expression operator/(const SymEngine::Expression expr, unsigned int i) {
    return expr / SymEngine::Expression(i);
}
inline SymEngine::Expression operator/(double d, const SymEngine::Expression expr) {
    return SymEngine::Expression(d) / expr;
}
inline SymEngine::Expression operator/(const SymEngine::Expression expr, double d) {
    return expr / SymEngine::Expression(d);
}

// --- operator* ---
inline SymEngine::Expression operator*(int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) * expr;
}
inline SymEngine::Expression operator*(const SymEngine::Expression expr, int i) {
    return expr * SymEngine::Expression(i);
}
inline SymEngine::Expression operator*(unsigned int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) * expr;
}
inline SymEngine::Expression operator*(const SymEngine::Expression expr, unsigned int i) {
    return expr * SymEngine::Expression(i);
}
inline SymEngine::Expression operator*(double d, const SymEngine::Expression expr) {
    return SymEngine::Expression(d) * expr;
}
inline SymEngine::Expression operator*(const SymEngine::Expression expr, double d) {
    return expr * SymEngine::Expression(d);
}

// --- operator+ ---
inline SymEngine::Expression operator+(int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) + expr;
}
inline SymEngine::Expression operator+(const SymEngine::Expression expr, int i) {
    return expr + SymEngine::Expression(i);
}
inline SymEngine::Expression operator+(unsigned int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) + expr;
}
inline SymEngine::Expression operator+(const SymEngine::Expression expr, unsigned int i) {
    return expr + SymEngine::Expression(i);
}
inline SymEngine::Expression operator+(double d, const SymEngine::Expression expr) {
    return SymEngine::Expression(d) + expr;
}
inline SymEngine::Expression operator+(const SymEngine::Expression expr, double d) {
    return expr + SymEngine::Expression(d);
}

// --- operator- ---
inline SymEngine::Expression operator-(int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) - expr;
}
inline SymEngine::Expression operator-(const SymEngine::Expression expr, int i) {
    return expr - SymEngine::Expression(i);
}
inline SymEngine::Expression operator-(unsigned int i, const SymEngine::Expression expr) {
    return SymEngine::Expression(i) - expr;
}
inline SymEngine::Expression operator-(const SymEngine::Expression expr, unsigned int i) {
    return expr - SymEngine::Expression(i);
}
inline SymEngine::Expression operator-(double d, const SymEngine::Expression expr) {
    return SymEngine::Expression(d) - expr;
}
inline SymEngine::Expression operator-(const SymEngine::Expression expr, double d) {
    return expr - SymEngine::Expression(d);
}

#endif // DA_WITH_SYMBOLIC
