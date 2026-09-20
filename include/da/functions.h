/**
 * @file functions.h
 * @brief Math functions for DAVector<double>.
 *
 * @details Ports of the math functions from ref/tpsa/src/da.cc:
 *   sqrt, exp, log, sin, cos, tan, asin, acos, atan,
 *   sinh, cosh, tanh, asinh, acosh, atanh,
 *   pow(int), pow(double), abs, erf, atan2.
 *
 * Complex DA variants (exp, sqrt, log, asin, acos, atan, asinh, acosh,
 * atanh, pow) are also provided.
 */

#pragma once

#include "da/davector.h"

#include <complex>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace da {

// ===========================================================================
// Real DA functions
// ===========================================================================

NDA sqrt(const NDA& v);
NDA exp(const NDA& v);
NDA log(const NDA& v);
NDA sin(const NDA& v);
NDA cos(const NDA& v);
NDA tan(const NDA& v);
NDA asin(const NDA& v);
NDA acos(const NDA& v);
NDA atan(const NDA& v);
NDA sinh(const NDA& v);
NDA cosh(const NDA& v);
NDA tanh(const NDA& v);
NDA asinh(const NDA& v);
NDA acosh(const NDA& v);
NDA atanh(const NDA& v);
NDA pow(const NDA& v, int order);
NDA pow(const NDA& v, double order);
double abs(const NDA& v);
NDA erf(const NDA& v);
NDA atan2(const NDA& y, const NDA& x);

// ===========================================================================
// Complex DA functions — template declarations over DAVector<T>
// Explicitly instantiated for T=double (always) and T=Expression (symbolic).
// ===========================================================================

template<class T> std::complex<DAVector<T>> exp(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> sqrt(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> log(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> asin(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> acos(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> atan(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> asinh(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> acosh(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> atanh(const std::complex<DAVector<T>>& c);
template<class T> std::complex<DAVector<T>> pow(const std::complex<DAVector<T>>& v, int order);
template<class T> std::complex<DAVector<T>> pow(const std::complex<DAVector<T>>& v, double order);

// Numerical abs: returns double (norm-like measure), unchanged.
double abs(const std::complex<NDA>& v);

#ifdef DA_WITH_SYMBOLIC

/// @brief Symbolic complex magnitude: sqrt(re^2 + im^2) as an SDA.
SDA abs(const std::complex<SDA>& v);

#endif // DA_WITH_SYMBOLIC (inner block for symbolic abs)

#ifdef DA_WITH_SYMBOLIC

// ===========================================================================
// Symbolic DA functions (SDA = DAVector<SymEngine::Expression>)
// Subset matching ref/tpsa_sym/include/sda.h
// ===========================================================================

SDA sqrt(const SDA& v);
SDA exp(const SDA& v);
SDA log(const SDA& v);
SDA sin(const SDA& v);
SDA cos(const SDA& v);
SDA tan(const SDA& v);
SDA asin(const SDA& v);
SDA acos(const SDA& v);
SDA atan(const SDA& v);
SDA sinh(const SDA& v);
SDA cosh(const SDA& v);
SDA tanh(const SDA& v);
SDA pow(const SDA& v, int order);
SDA pow(const SDA& v, double order);
SDA erf(const SDA& v);

#endif // DA_WITH_SYMBOLIC

} // namespace da
