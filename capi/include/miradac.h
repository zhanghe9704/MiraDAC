/* miradac.h — C API of the MiraDAC differential algebra library.
 *
 * Conventions (DEVELOPMENT_PLAN_JULIA.md A.3):
 * - Functions that can fail return an mdac_status and pass results through
 *   out-parameters; mdac_last_error() then gives the message.
 * - Every DA vector handle is owned by the caller and released with its
 *   mdac_*_free function, which never fails (also after its env was closed).
 * - Env handles are not owned: envs belong to the library.
 * - Variable-length outputs follow the size-query convention: the function
 *   writes at most cap items to buf and always sets *n to the full size.
 * - Indices are 0-based.
 */
#ifndef MIRADAC_H
#define MIRADAC_H

#include <stddef.h>
#include <stdint.h>

#define MDAC_API __attribute__((visibility("default")))

/* Bumped whenever the C API changes incompatibly. */
#define MDAC_ABI_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

typedef int mdac_status;
enum {
    MDAC_OK = 0,
    MDAC_ERR_ENV = 1,         /* env mismatch, cleared env, no env */
    MDAC_ERR_VALUE = 2,       /* invalid argument, math domain error */
    MDAC_ERR_INDEX = 3,       /* index out of range */
    MDAC_ERR_POOL = 4,        /* the env's pool has no free vector */
    MDAC_ERR_RUNTIME = 5,     /* anything else */
    MDAC_ERR_UNSUPPORTED = 6  /* a symbolic function in a numeric-only build */
};

typedef struct mdac_env mdac_env;
typedef struct mdac_nda mdac_nda;

MDAC_API int mdac_abi_version(void);
MDAC_API const char* mdac_version(void);

/* Message of the last failure on the calling thread. */
MDAC_API const char* mdac_last_error(void);

/* ---- Envs ---------------------------------------------------------------- */

/* The default env (C++ da_init / da_clear); init also makes it current. */
MDAC_API mdac_status mdac_init(unsigned order, unsigned nvars, unsigned poolsize, int table);
MDAC_API mdac_status mdac_clear(void);

MDAC_API mdac_status mdac_env_current(mdac_env** out);
MDAC_API mdac_status mdac_env_default(mdac_env** out);
/* A new env; the current env does not change. */
MDAC_API mdac_status mdac_env_make(unsigned order, unsigned nvars, unsigned poolsize, int table,
                                   mdac_env** out);
MDAC_API mdac_status mdac_env_select(mdac_env* e);
/* Makes e (may be null) current and returns the previous current env. */
MDAC_API mdac_status mdac_env_exchange(mdac_env* e, mdac_env** prev);
/* Releases the pool; the env's vectors then fail with MDAC_ERR_ENV. Idempotent. */
MDAC_API mdac_status mdac_env_close(mdac_env* e);

/* order: the current truncation order; max_order: the order at creation. */
MDAC_API mdac_status mdac_env_order(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_max_order(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_nvars(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_full_length(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_poolsize(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_count(mdac_env* e, unsigned* out);
MDAC_API mdac_status mdac_env_remain(mdac_env* e, unsigned* out);
MDAC_API int mdac_env_retired(const mdac_env* e);
/* *ok = 0 if n is larger than max_order. */
MDAC_API mdac_status mdac_env_change_order(mdac_env* e, unsigned n, int* ok);
MDAC_API mdac_status mdac_env_restore_order(mdac_env* e);

/* Coefficients below eps in magnitude are dropped (shared by all envs). */
MDAC_API double mdac_get_eps(void);
MDAC_API mdac_status mdac_set_eps(double x);

/* ---- NDA: lifecycle and inspection --------------------------------------- */

MDAC_API mdac_status mdac_nda_new(mdac_env* env, double x, mdac_nda** out);
MDAC_API mdac_status mdac_nda_from_coeffs(mdac_env* env, const double* c, size_t n, mdac_nda** out);
MDAC_API mdac_status mdac_nda_var(mdac_env* env, unsigned i, mdac_nda** out);
MDAC_API mdac_status mdac_nda_copy(const mdac_nda* v, mdac_nda** out);
MDAC_API void mdac_nda_free(mdac_nda* v);
MDAC_API mdac_env* mdac_nda_env(const mdac_nda* v);
/* A copy of v in env (C++ da::import_to); the current env does not change.
 * MDAC_ERR_VALUE if the two envs have different variables or orders. The
 * other types' _import functions work the same way. */
MDAC_API mdac_status mdac_nda_import(mdac_env* env, const mdac_nda* v, mdac_nda** out);

MDAC_API mdac_status mdac_nda_con(const mdac_nda* v, double* out);
/* Resets v to the constant x. */
MDAC_API mdac_status mdac_nda_set_con(mdac_nda* v, double x);
/* Number of stored coefficients (up to the last non-zero one). */
MDAC_API mdac_status mdac_nda_length(const mdac_nda* v, size_t* out);
/* Number of non-zero coefficients. */
MDAC_API mdac_status mdac_nda_nterms(const mdac_nda* v, size_t* out);
/* The largest coefficient magnitude. */
MDAC_API mdac_status mdac_nda_norm(const mdac_nda* v, double* out);
/* The stored coefficients in monomial order (size-query convention). */
MDAC_API mdac_status mdac_nda_coeffs(const mdac_nda* v, double* buf, size_t cap, size_t* n);
/* The coefficient of the monomial with exponents exps[0..k). */
MDAC_API mdac_status mdac_nda_coeff(const mdac_nda* v, const int* exps, size_t k, double* out);
MDAC_API mdac_status mdac_nda_set_coeff(mdac_nda* v, const int* exps, size_t k, double x);
/* Exponents (nvars entries) and value of coefficient i in monomial order. */
MDAC_API mdac_status mdac_nda_index_term(const mdac_nda* v, size_t i, int* exps, double* out);
MDAC_API mdac_status mdac_nda_iszero(const mdac_nda* v, double eps, int* out);
MDAC_API mdac_status mdac_nda_clean(mdac_nda* v, double eps);
MDAC_API mdac_status mdac_nda_reset(mdac_nda* v);
/* The C++ operator<< text; *n counts the terminating NUL, which is always
 * written when cap > 0 (size-query convention). */
MDAC_API mdac_status mdac_nda_to_string(const mdac_nda* v, char* buf, size_t cap, size_t* n);

/* ---- NDA: arithmetic ------------------------------------------------------
 * mdac_nda_<op> allocates the result; mdac_nda_<op>_into writes it into the
 * existing vector out (same env; out may be an operand). _d: a double on the
 * right; d<op>: a double on the left. */

#define MDAC_NDA_BINOP_DECL(op)                                                                   \
    MDAC_API mdac_status mdac_nda_##op(const mdac_nda* a, const mdac_nda* b, mdac_nda** out);     \
    MDAC_API mdac_status mdac_nda_##op##_d(const mdac_nda* a, double x, mdac_nda** out);          \
    MDAC_API mdac_status mdac_nda_d##op(double x, const mdac_nda* a, mdac_nda** out);             \
    MDAC_API mdac_status mdac_nda_##op##_into(mdac_nda* out, const mdac_nda* a, const mdac_nda* b); \
    MDAC_API mdac_status mdac_nda_##op##_d_into(mdac_nda* out, const mdac_nda* a, double x);      \
    MDAC_API mdac_status mdac_nda_d##op##_into(mdac_nda* out, double x, const mdac_nda* a);

MDAC_NDA_BINOP_DECL(add)
MDAC_NDA_BINOP_DECL(sub)
MDAC_NDA_BINOP_DECL(mul)
MDAC_NDA_BINOP_DECL(div)

MDAC_API mdac_status mdac_nda_neg(const mdac_nda* a, mdac_nda** out);
MDAC_API mdac_status mdac_nda_neg_into(mdac_nda* out, const mdac_nda* a);
MDAC_API mdac_status mdac_nda_pow_i(const mdac_nda* a, int n, mdac_nda** out);
MDAC_API mdac_status mdac_nda_pow_i_into(mdac_nda* out, const mdac_nda* a, int n);
MDAC_API mdac_status mdac_nda_pow_d(const mdac_nda* a, double x, mdac_nda** out);
MDAC_API mdac_status mdac_nda_pow_d_into(mdac_nda* out, const mdac_nda* a, double x);

/* ---- NDA: math functions (mdac_nda_exp(a, out), mdac_nda_exp_into(out, a)) - */

#define MDAC_NDA_FUNC_DECL(f)                                                  \
    MDAC_API mdac_status mdac_nda_##f(const mdac_nda* a, mdac_nda** out);      \
    MDAC_API mdac_status mdac_nda_##f##_into(mdac_nda* out, const mdac_nda* a);

MDAC_NDA_FUNC_DECL(sqrt)
MDAC_NDA_FUNC_DECL(exp)
MDAC_NDA_FUNC_DECL(log)
MDAC_NDA_FUNC_DECL(sin)
MDAC_NDA_FUNC_DECL(cos)
MDAC_NDA_FUNC_DECL(tan)
MDAC_NDA_FUNC_DECL(asin)
MDAC_NDA_FUNC_DECL(acos)
MDAC_NDA_FUNC_DECL(atan)
MDAC_NDA_FUNC_DECL(sinh)
MDAC_NDA_FUNC_DECL(cosh)
MDAC_NDA_FUNC_DECL(tanh)
MDAC_NDA_FUNC_DECL(asinh)
MDAC_NDA_FUNC_DECL(acosh)
MDAC_NDA_FUNC_DECL(atanh)
MDAC_NDA_FUNC_DECL(erf)

/* abs(NDA): the norm. */
MDAC_API mdac_status mdac_nda_abs(const mdac_nda* a, double* out);

/* ---- NDA lists --------------------------------------------------------------
 * A list holds copies of vectors of one env (MDAC_ERR_ENV otherwise). */

typedef struct mdac_ndalist mdac_ndalist;

MDAC_API mdac_status mdac_ndalist_new(mdac_ndalist** out);
/* A list of copies of vs[0..n). */
MDAC_API mdac_status mdac_ndalist_from(const mdac_nda* const* vs, size_t n, mdac_ndalist** out);
MDAC_API void mdac_ndalist_free(mdac_ndalist* l);
MDAC_API size_t mdac_ndalist_length(const mdac_ndalist* l);
/* The env of the list's vectors; null for an empty list. */
MDAC_API mdac_env* mdac_ndalist_env(const mdac_ndalist* l);
/* A copy of element i. */
MDAC_API mdac_status mdac_ndalist_get(const mdac_ndalist* l, size_t i, mdac_nda** out);
/* Copies v into element i. */
MDAC_API mdac_status mdac_ndalist_set(mdac_ndalist* l, size_t i, const mdac_nda* v);
/* Appends a copy of v. */
MDAC_API mdac_status mdac_ndalist_push(mdac_ndalist* l, const mdac_nda* v);

/* ---- NDA algorithms ---------------------------------------------------------
 * Variables are 0-based; an index >= nvars gives MDAC_ERR_INDEX. */

MDAC_API mdac_status mdac_nda_der(const mdac_nda* v, unsigned i, mdac_nda** out);
MDAC_API mdac_status mdac_nda_integ(const mdac_nda* v, unsigned i, mdac_nda** out);
/* iv with x (a number, a vector) substituted for variable i. */
MDAC_API mdac_status mdac_nda_substitute_d(const mdac_nda* iv, unsigned i, double x, mdac_nda** out);
MDAC_API mdac_status mdac_nda_substitute(const mdac_nda* iv, unsigned i, const mdac_nda* v,
                                         mdac_nda** out);
/* v[j] substituted for variable ids[j], j < k (k == length of v, distinct ids). */
MDAC_API mdac_status mdac_nda_substitute_multi(const mdac_nda* iv, const unsigned* ids, size_t k,
                                               const mdac_ndalist* v, mdac_nda** out);
MDAC_API mdac_status mdac_ndalist_substitute(const mdac_ndalist* m, const unsigned* ids, size_t k,
                                             const mdac_ndalist* v, mdac_ndalist** out);
/* m[i](v[0], ..., v[nvars-1]); v has nvars vectors. */
MDAC_API mdac_status mdac_ndalist_compose(const mdac_ndalist* m, const mdac_ndalist* v,
                                          mdac_ndalist** out);
/* The values of m at the point pt (n == nvars); out has length(m) entries. */
MDAC_API mdac_status mdac_ndalist_compose_d(const mdac_ndalist* m, const double* pt, size_t n,
                                            double* out);
/* The same at a complex point: pt holds n (re, im) pairs, out length(m) pairs. */
MDAC_API mdac_status mdac_ndalist_compose_z(const mdac_ndalist* m, const double* pt, size_t n,
                                            double* out);
/* The inverse of the map m[0..dim) (zero constant parts), a new list of dim
 * vectors; MDAC_ERR_VALUE if the linear part is singular. */
MDAC_API mdac_status mdac_ndalist_inv_map(const mdac_ndalist* m, int dim, mdac_ndalist** out);
/* m at npts points stored one after another (nvars values each):
 * out[p * length(m) + c] = m[c] at point p, equal to mdac_ndalist_compose_d. */
MDAC_API mdac_status mdac_ndalist_evaluate_map(const mdac_ndalist* m, const double* pts,
                                               size_t npts, double* out);
/* v at the point pt (n == nvars). */
MDAC_API mdac_status mdac_nda_eval(const mdac_nda* v, const double* pt, size_t n, double* out);
/* The exponents of every monomial of env, nvars per monomial, in monomial
 * order (size-query convention; *n = full_length * nvars). */
MDAC_API mdac_status mdac_exponents(mdac_env* env, int* buf, size_t cap, size_t* n);

/* ---- CNDA: complex numeric DA vectors (C++ std::complex<NDA>) ----------------
 * A CNDA's real and imaginary parts belong to one env. */

typedef struct mdac_cnda mdac_cnda;

/* re + i*im from vectors of one env; im may be null (a zero part). */
MDAC_API mdac_status mdac_cnda_new(const mdac_nda* re, const mdac_nda* im, mdac_cnda** out);
/* The constant re + i*im. */
MDAC_API mdac_status mdac_cnda_new_z(mdac_env* env, double re, double im, mdac_cnda** out);
MDAC_API mdac_status mdac_cnda_copy(const mdac_cnda* v, mdac_cnda** out);
MDAC_API void mdac_cnda_free(mdac_cnda* v);
MDAC_API mdac_env* mdac_cnda_env(const mdac_cnda* v);
MDAC_API mdac_status mdac_cnda_import(mdac_env* env, const mdac_cnda* v, mdac_cnda** out);
/* Copies of the parts. */
MDAC_API mdac_status mdac_cnda_real(const mdac_cnda* v, mdac_nda** out);
MDAC_API mdac_status mdac_cnda_imag(const mdac_cnda* v, mdac_nda** out);
/* Copy x into a part. */
MDAC_API mdac_status mdac_cnda_set_real(mdac_cnda* v, const mdac_nda* x);
MDAC_API mdac_status mdac_cnda_set_imag(mdac_cnda* v, const mdac_nda* x);
/* The C++ operator<< text, as mdac_nda_to_string. */
MDAC_API mdac_status mdac_cnda_to_string(const mdac_cnda* v, char* buf, size_t cap, size_t* n);

/* ---- CNDA: arithmetic -------------------------------------------------------
 * Operands: c a CNDA, n an NDA, d a double, z a complex (re, im). The suffix
 * names the right operand, the prefix the left one: mdac_cnda_<op> (c, c),
 * _n/n<op> (c, n)/(n, c), _d/d<op>, _z/z<op>; mdac_nda_<op>_z and
 * mdac_nda_z<op> (n, z)/(z, n) give a CNDA too. _into forms write into the
 * existing CNDA out (same env; out may be an operand). */

#define MDAC_CNDA_BINOP_DECL(op)                                                                        \
    MDAC_API mdac_status mdac_cnda_##op(const mdac_cnda* a, const mdac_cnda* b, mdac_cnda** out);       \
    MDAC_API mdac_status mdac_cnda_##op##_n(const mdac_cnda* a, const mdac_nda* b, mdac_cnda** out);    \
    MDAC_API mdac_status mdac_cnda_n##op(const mdac_nda* a, const mdac_cnda* b, mdac_cnda** out);       \
    MDAC_API mdac_status mdac_cnda_##op##_d(const mdac_cnda* a, double x, mdac_cnda** out);             \
    MDAC_API mdac_status mdac_cnda_d##op(double x, const mdac_cnda* a, mdac_cnda** out);                \
    MDAC_API mdac_status mdac_cnda_##op##_z(const mdac_cnda* a, double re, double im, mdac_cnda** out); \
    MDAC_API mdac_status mdac_cnda_z##op(double re, double im, const mdac_cnda* a, mdac_cnda** out);    \
    MDAC_API mdac_status mdac_nda_##op##_z(const mdac_nda* a, double re, double im, mdac_cnda** out);   \
    MDAC_API mdac_status mdac_nda_z##op(double re, double im, const mdac_nda* a, mdac_cnda** out);      \
    MDAC_API mdac_status mdac_cnda_##op##_into(mdac_cnda* out, const mdac_cnda* a, const mdac_cnda* b); \
    MDAC_API mdac_status mdac_cnda_##op##_n_into(mdac_cnda* out, const mdac_cnda* a, const mdac_nda* b); \
    MDAC_API mdac_status mdac_cnda_n##op##_into(mdac_cnda* out, const mdac_nda* a, const mdac_cnda* b); \
    MDAC_API mdac_status mdac_cnda_##op##_d_into(mdac_cnda* out, const mdac_cnda* a, double x);         \
    MDAC_API mdac_status mdac_cnda_d##op##_into(mdac_cnda* out, double x, const mdac_cnda* a);          \
    MDAC_API mdac_status mdac_cnda_##op##_z_into(mdac_cnda* out, const mdac_cnda* a, double re, double im); \
    MDAC_API mdac_status mdac_cnda_z##op##_into(mdac_cnda* out, double re, double im, const mdac_cnda* a); \
    MDAC_API mdac_status mdac_nda_##op##_z_into(mdac_cnda* out, const mdac_nda* a, double re, double im); \
    MDAC_API mdac_status mdac_nda_z##op##_into(mdac_cnda* out, double re, double im, const mdac_nda* a);

MDAC_CNDA_BINOP_DECL(add)
MDAC_CNDA_BINOP_DECL(sub)
MDAC_CNDA_BINOP_DECL(mul)
MDAC_CNDA_BINOP_DECL(div)

MDAC_API mdac_status mdac_cnda_neg(const mdac_cnda* a, mdac_cnda** out);
MDAC_API mdac_status mdac_cnda_neg_into(mdac_cnda* out, const mdac_cnda* a);
MDAC_API mdac_status mdac_cnda_pow_i(const mdac_cnda* a, int n, mdac_cnda** out);
MDAC_API mdac_status mdac_cnda_pow_i_into(mdac_cnda* out, const mdac_cnda* a, int n);
MDAC_API mdac_status mdac_cnda_pow_d(const mdac_cnda* a, double x, mdac_cnda** out);
MDAC_API mdac_status mdac_cnda_pow_d_into(mdac_cnda* out, const mdac_cnda* a, double x);

/* ---- CNDA: math functions (C++ has no complex sin cos tan sinh cosh tanh erf) */

#define MDAC_CNDA_FUNC_DECL(f)                                                   \
    MDAC_API mdac_status mdac_cnda_##f(const mdac_cnda* a, mdac_cnda** out);     \
    MDAC_API mdac_status mdac_cnda_##f##_into(mdac_cnda* out, const mdac_cnda* a);

MDAC_CNDA_FUNC_DECL(sqrt)
MDAC_CNDA_FUNC_DECL(exp)
MDAC_CNDA_FUNC_DECL(log)
MDAC_CNDA_FUNC_DECL(asin)
MDAC_CNDA_FUNC_DECL(acos)
MDAC_CNDA_FUNC_DECL(atan)
MDAC_CNDA_FUNC_DECL(asinh)
MDAC_CNDA_FUNC_DECL(acosh)
MDAC_CNDA_FUNC_DECL(atanh)

/* abs(CNDA): the larger of the parts' norms. */
MDAC_API mdac_status mdac_cnda_abs(const mdac_cnda* a, double* out);

/* ---- CNDA lists and composition ----------------------------------------------
 * As NDA lists. The compose functions are C++ cd_composition, giving a new
 * list of length(m): m[i](v[0], ..., v[nvars-1]) with v of nvars vectors.
 * An NDA map at complex points is mdac_ndalist_compose_z. */

typedef struct mdac_cndalist mdac_cndalist;

MDAC_API mdac_status mdac_cndalist_new(mdac_cndalist** out);
MDAC_API mdac_status mdac_cndalist_from(const mdac_cnda* const* vs, size_t n, mdac_cndalist** out);
MDAC_API void mdac_cndalist_free(mdac_cndalist* l);
MDAC_API size_t mdac_cndalist_length(const mdac_cndalist* l);
MDAC_API mdac_env* mdac_cndalist_env(const mdac_cndalist* l);
MDAC_API mdac_status mdac_cndalist_get(const mdac_cndalist* l, size_t i, mdac_cnda** out);
MDAC_API mdac_status mdac_cndalist_set(mdac_cndalist* l, size_t i, const mdac_cnda* v);
MDAC_API mdac_status mdac_cndalist_push(mdac_cndalist* l, const mdac_cnda* v);

/* An NDA map at CNDA arguments. */
MDAC_API mdac_status mdac_ndalist_compose_c(const mdac_ndalist* m, const mdac_cndalist* v,
                                            mdac_cndalist** out);
/* A CNDA map at CNDA arguments. */
MDAC_API mdac_status mdac_cndalist_compose(const mdac_cndalist* m, const mdac_cndalist* v,
                                           mdac_cndalist** out);
/* A CNDA map at NDA arguments. */
MDAC_API mdac_status mdac_cndalist_compose_n(const mdac_cndalist* m, const mdac_ndalist* v,
                                             mdac_cndalist** out);

/* ---- Symbolic support ---------------------------------------------------------
 * Every function below exists in every build. In a build without symbolic
 * support (WITH_SYMBOLIC=OFF) those returning a status give
 * MDAC_ERR_UNSUPPORTED, the free functions do nothing, and the length and env
 * functions return 0 and null. */

/* 1 if the library was built with symbolic support, else 0. */
MDAC_API int mdac_has_symbolic(void);

/* ---- Expr: symbolic scalars (C++ SymEngine::Expression) ----------------------
 * Not tied to an env. Results are new handles, freed with mdac_expr_free. */

typedef struct mdac_expr mdac_expr;

MDAC_API mdac_status mdac_expr_new_d(double x, mdac_expr** out);
/* An exact integer. */
MDAC_API mdac_status mdac_expr_new_i(int64_t x, mdac_expr** out);
/* SymEngine::parse(s), e.g. "a + 2*b"; MDAC_ERR_VALUE for bad input. */
MDAC_API mdac_status mdac_expr_parse(const char* s, mdac_expr** out);
MDAC_API mdac_status mdac_expr_symbol(const char* name, mdac_expr** out);
MDAC_API mdac_status mdac_expr_copy(const mdac_expr* x, mdac_expr** out);
MDAC_API void mdac_expr_free(mdac_expr* x);
/* The text of x, as mdac_nda_to_string. */
MDAC_API mdac_status mdac_expr_to_string(const mdac_expr* x, char* buf, size_t cap, size_t* n);
/* The numeric value; MDAC_ERR_VALUE if free symbols remain. */
MDAC_API mdac_status mdac_expr_to_double(const mdac_expr* x, double* out);

/* mdac_expr_<op>(a, b), _d: a double on the right, d<op>: on the left. pow is a**b. */
#define MDAC_EXPR_BINOP_DECL(op)                                                                \
    MDAC_API mdac_status mdac_expr_##op(const mdac_expr* a, const mdac_expr* b, mdac_expr** out); \
    MDAC_API mdac_status mdac_expr_##op##_d(const mdac_expr* a, double x, mdac_expr** out);      \
    MDAC_API mdac_status mdac_expr_d##op(double x, const mdac_expr* a, mdac_expr** out);

MDAC_EXPR_BINOP_DECL(add)
MDAC_EXPR_BINOP_DECL(sub)
MDAC_EXPR_BINOP_DECL(mul)
MDAC_EXPR_BINOP_DECL(div)
MDAC_EXPR_BINOP_DECL(pow)

MDAC_API mdac_status mdac_expr_neg(const mdac_expr* a, mdac_expr** out);
/* Structural equality, and SymEngine's hash, consistent with it. */
MDAC_API mdac_status mdac_expr_eq(const mdac_expr* a, const mdac_expr* b, int* out);
MDAC_API mdac_status mdac_expr_hash(const mdac_expr* a, uint64_t* out);
/* 1 if x is zero (the test SDA uses: simplify and expand). */
MDAC_API mdac_status mdac_expr_is_zero(const mdac_expr* x, int* out);
/* x with vals[j] substituted for keys[j], j < n. */
MDAC_API mdac_status mdac_expr_subs(const mdac_expr* x, size_t n, const mdac_expr* const* keys,
                                    const mdac_expr* const* vals, mdac_expr** out);
MDAC_API mdac_status mdac_expr_expand(const mdac_expr* x, mdac_expr** out);
/* The derivative with respect to the symbol sym (MDAC_ERR_VALUE if not a symbol). */
MDAC_API mdac_status mdac_expr_diff(const mdac_expr* x, const mdac_expr* sym, mdac_expr** out);
/* The simplification SDA uses (mdac_sda_simplify). */
MDAC_API mdac_status mdac_expr_simplify(const mdac_expr* x, mdac_expr** out);
/* The symbols of x as new handles (size-query convention). */
MDAC_API mdac_status mdac_expr_free_symbols(const mdac_expr* x, mdac_expr** buf, size_t cap,
                                            size_t* n);

/* ---- SDA: symbolic DA vectors (C++ DAVector<SymEngine::Expression>) -----------
 * As NDA. Coefficients come and go as new mdac_expr handles. Not provided for
 * SDA (C++ has no SDA version, or the kernel ignores symbolic coefficients):
 * from_coeffs, norm, abs, the eps argument of iszero and clean, eval at a
 * point, asinh, acosh, atanh, inv_map, evaluate_map. */

typedef struct mdac_sda mdac_sda;
typedef struct mdac_sdalist mdac_sdalist;

/* The constant x (null: zero). */
MDAC_API mdac_status mdac_sda_new(mdac_env* env, const mdac_expr* x, mdac_sda** out);
MDAC_API mdac_status mdac_sda_new_d(mdac_env* env, double x, mdac_sda** out);
/* The symbolic base vector of variable i (promote(var(i))). */
MDAC_API mdac_status mdac_sda_var(mdac_env* env, unsigned i, mdac_sda** out);
/* v as an SDA of its env; integer-valued coefficients become exact integers. */
MDAC_API mdac_status mdac_sda_promote(const mdac_nda* v, mdac_sda** out);
MDAC_API mdac_status mdac_sda_copy(const mdac_sda* v, mdac_sda** out);
MDAC_API void mdac_sda_free(mdac_sda* v);
MDAC_API mdac_env* mdac_sda_env(const mdac_sda* v);
MDAC_API mdac_status mdac_sda_import(mdac_env* env, const mdac_sda* v, mdac_sda** out);
/* mdac_sda_promote(v) with the result in env (C++ da::promote_to); the
 * current env does not change. MDAC_ERR_VALUE if the layouts differ. */
MDAC_API mdac_status mdac_nda_promote_to(mdac_env* env, const mdac_nda* v, mdac_sda** out);

MDAC_API mdac_status mdac_sda_con(const mdac_sda* v, mdac_expr** out);
/* Resets v to the constant x. */
MDAC_API mdac_status mdac_sda_set_con(mdac_sda* v, const mdac_expr* x);
MDAC_API mdac_status mdac_sda_length(const mdac_sda* v, size_t* out);
MDAC_API mdac_status mdac_sda_nterms(const mdac_sda* v, size_t* out);
/* The stored coefficients as new handles (size-query convention). */
MDAC_API mdac_status mdac_sda_coeffs(const mdac_sda* v, mdac_expr** buf, size_t cap, size_t* n);
MDAC_API mdac_status mdac_sda_coeff(const mdac_sda* v, const int* exps, size_t k, mdac_expr** out);
MDAC_API mdac_status mdac_sda_set_coeff(mdac_sda* v, const int* exps, size_t k, const mdac_expr* x);
MDAC_API mdac_status mdac_sda_index_term(const mdac_sda* v, size_t i, int* exps, mdac_expr** out);
/* Zero tests use mdac_expr_is_zero. */
MDAC_API mdac_status mdac_sda_iszero(const mdac_sda* v, int* out);
MDAC_API mdac_status mdac_sda_clean(mdac_sda* v);
MDAC_API mdac_status mdac_sda_reset(mdac_sda* v);
MDAC_API mdac_status mdac_sda_to_string(const mdac_sda* v, char* buf, size_t cap, size_t* n);

/* ---- SDA: arithmetic ----------------------------------------------------------
 * Operands: s an SDA, n an NDA, e an Expr, d a double; the suffix names the
 * right operand, the prefix the left one, as for CNDA: mdac_sda_<op> (s, s),
 * _n/n<op>, _e/e<op>, _d/d<op>; mdac_nda_<op>_e and mdac_nda_e<op> (n, e)/(e, n)
 * give an SDA too. _into forms write into the existing SDA out (same env;
 * out may be an operand). */

#define MDAC_SDA_BINOP_DECL(op)                                                                       \
    MDAC_API mdac_status mdac_sda_##op(const mdac_sda* a, const mdac_sda* b, mdac_sda** out);         \
    MDAC_API mdac_status mdac_sda_##op##_n(const mdac_sda* a, const mdac_nda* b, mdac_sda** out);     \
    MDAC_API mdac_status mdac_sda_n##op(const mdac_nda* a, const mdac_sda* b, mdac_sda** out);        \
    MDAC_API mdac_status mdac_sda_##op##_e(const mdac_sda* a, const mdac_expr* x, mdac_sda** out);    \
    MDAC_API mdac_status mdac_sda_e##op(const mdac_expr* x, const mdac_sda* a, mdac_sda** out);       \
    MDAC_API mdac_status mdac_sda_##op##_d(const mdac_sda* a, double x, mdac_sda** out);              \
    MDAC_API mdac_status mdac_sda_d##op(double x, const mdac_sda* a, mdac_sda** out);                 \
    MDAC_API mdac_status mdac_nda_##op##_e(const mdac_nda* a, const mdac_expr* x, mdac_sda** out);    \
    MDAC_API mdac_status mdac_nda_e##op(const mdac_expr* x, const mdac_nda* a, mdac_sda** out);       \
    MDAC_API mdac_status mdac_sda_##op##_into(mdac_sda* out, const mdac_sda* a, const mdac_sda* b);   \
    MDAC_API mdac_status mdac_sda_##op##_n_into(mdac_sda* out, const mdac_sda* a, const mdac_nda* b); \
    MDAC_API mdac_status mdac_sda_n##op##_into(mdac_sda* out, const mdac_nda* a, const mdac_sda* b);  \
    MDAC_API mdac_status mdac_sda_##op##_e_into(mdac_sda* out, const mdac_sda* a, const mdac_expr* x); \
    MDAC_API mdac_status mdac_sda_e##op##_into(mdac_sda* out, const mdac_expr* x, const mdac_sda* a); \
    MDAC_API mdac_status mdac_sda_##op##_d_into(mdac_sda* out, const mdac_sda* a, double x);         \
    MDAC_API mdac_status mdac_sda_d##op##_into(mdac_sda* out, double x, const mdac_sda* a);          \
    MDAC_API mdac_status mdac_nda_##op##_e_into(mdac_sda* out, const mdac_nda* a, const mdac_expr* x); \
    MDAC_API mdac_status mdac_nda_e##op##_into(mdac_sda* out, const mdac_expr* x, const mdac_nda* a);

MDAC_SDA_BINOP_DECL(add)
MDAC_SDA_BINOP_DECL(sub)
MDAC_SDA_BINOP_DECL(mul)
MDAC_SDA_BINOP_DECL(div)

MDAC_API mdac_status mdac_sda_neg(const mdac_sda* a, mdac_sda** out);
MDAC_API mdac_status mdac_sda_neg_into(mdac_sda* out, const mdac_sda* a);
MDAC_API mdac_status mdac_sda_pow_i(const mdac_sda* a, int n, mdac_sda** out);
MDAC_API mdac_status mdac_sda_pow_i_into(mdac_sda* out, const mdac_sda* a, int n);
MDAC_API mdac_status mdac_sda_pow_d(const mdac_sda* a, double x, mdac_sda** out);
MDAC_API mdac_status mdac_sda_pow_d_into(mdac_sda* out, const mdac_sda* a, double x);

/* ---- SDA: math functions (C++ has no SDA asinh acosh atanh) ------------------ */

#define MDAC_SDA_FUNC_DECL(f)                                                  \
    MDAC_API mdac_status mdac_sda_##f(const mdac_sda* a, mdac_sda** out);      \
    MDAC_API mdac_status mdac_sda_##f##_into(mdac_sda* out, const mdac_sda* a);

MDAC_SDA_FUNC_DECL(sqrt)
MDAC_SDA_FUNC_DECL(exp)
MDAC_SDA_FUNC_DECL(log)
MDAC_SDA_FUNC_DECL(sin)
MDAC_SDA_FUNC_DECL(cos)
MDAC_SDA_FUNC_DECL(tan)
MDAC_SDA_FUNC_DECL(asin)
MDAC_SDA_FUNC_DECL(acos)
MDAC_SDA_FUNC_DECL(atan)
MDAC_SDA_FUNC_DECL(sinh)
MDAC_SDA_FUNC_DECL(cosh)
MDAC_SDA_FUNC_DECL(tanh)
MDAC_SDA_FUNC_DECL(erf)

/* ---- SDA: per-coefficient operations and evaluation -------------------------- */

/* A new SDA with every coefficient simplified, expanded, or with vals[j]
 * substituted for keys[j] (j < n). */
MDAC_API mdac_status mdac_sda_simplify(const mdac_sda* v, mdac_sda** out);
MDAC_API mdac_status mdac_sda_expand(const mdac_sda* v, mdac_sda** out);
MDAC_API mdac_status mdac_sda_subs(const mdac_sda* v, size_t n, const mdac_expr* const* keys,
                                   const mdac_expr* const* vals, mdac_sda** out);
/* The NDA of v with vals[j] for the symbol keys[j] (j < n); every symbol of v
 * needs a value (MDAC_ERR_VALUE otherwise). */
MDAC_API mdac_status mdac_sda_evaluate(const mdac_sda* v, size_t n, const mdac_expr* const* keys,
                                       const double* vals, mdac_nda** out);

/* ---- SDA lists and algorithms -------------------------------------------------
 * As NDA lists and algorithms. */

MDAC_API mdac_status mdac_sdalist_new(mdac_sdalist** out);
MDAC_API mdac_status mdac_sdalist_from(const mdac_sda* const* vs, size_t n, mdac_sdalist** out);
MDAC_API void mdac_sdalist_free(mdac_sdalist* l);
MDAC_API size_t mdac_sdalist_length(const mdac_sdalist* l);
MDAC_API mdac_env* mdac_sdalist_env(const mdac_sdalist* l);
MDAC_API mdac_status mdac_sdalist_get(const mdac_sdalist* l, size_t i, mdac_sda** out);
MDAC_API mdac_status mdac_sdalist_set(mdac_sdalist* l, size_t i, const mdac_sda* v);
MDAC_API mdac_status mdac_sdalist_push(mdac_sdalist* l, const mdac_sda* v);

MDAC_API mdac_status mdac_sda_der(const mdac_sda* v, unsigned i, mdac_sda** out);
MDAC_API mdac_status mdac_sda_integ(const mdac_sda* v, unsigned i, mdac_sda** out);
MDAC_API mdac_status mdac_sda_substitute_d(const mdac_sda* iv, unsigned i, double x, mdac_sda** out);
MDAC_API mdac_status mdac_sda_substitute(const mdac_sda* iv, unsigned i, const mdac_sda* v,
                                         mdac_sda** out);
MDAC_API mdac_status mdac_sda_substitute_multi(const mdac_sda* iv, const unsigned* ids, size_t k,
                                               const mdac_sdalist* v, mdac_sda** out);
MDAC_API mdac_status mdac_sdalist_substitute(const mdac_sdalist* m, const unsigned* ids, size_t k,
                                             const mdac_sdalist* v, mdac_sdalist** out);
MDAC_API mdac_status mdac_sdalist_compose(const mdac_sdalist* m, const mdac_sdalist* v,
                                          mdac_sdalist** out);
/* m at the point pt (n == nvars): length(m) new Expr handles in out. */
MDAC_API mdac_status mdac_sdalist_compose_d(const mdac_sdalist* m, const double* pt, size_t n,
                                            mdac_expr** out);

/* ---- CSDA: complex symbolic DA vectors (C++ std::complex<SDA>) ---------------
 * As CNDA, with SDA parts of one env. */

typedef struct mdac_csda mdac_csda;
typedef struct mdac_csdalist mdac_csdalist;

/* re + i*im (im null: zero). */
MDAC_API mdac_status mdac_csda_new(const mdac_sda* re, const mdac_sda* im, mdac_csda** out);
/* The constant re + i*im in env. */
MDAC_API mdac_status mdac_csda_new_z(mdac_env* env, double re, double im, mdac_csda** out);
/* v as a CSDA of its env (both parts promoted as by mdac_sda_promote). */
MDAC_API mdac_status mdac_csda_promote(const mdac_cnda* v, mdac_csda** out);
MDAC_API mdac_status mdac_csda_copy(const mdac_csda* v, mdac_csda** out);
MDAC_API void mdac_csda_free(mdac_csda* v);
MDAC_API mdac_env* mdac_csda_env(const mdac_csda* v);
MDAC_API mdac_status mdac_csda_import(mdac_env* env, const mdac_csda* v, mdac_csda** out);
/* Copies of the parts. */
MDAC_API mdac_status mdac_csda_real(const mdac_csda* v, mdac_sda** out);
MDAC_API mdac_status mdac_csda_imag(const mdac_csda* v, mdac_sda** out);
MDAC_API mdac_status mdac_csda_set_real(mdac_csda* v, const mdac_sda* x);
MDAC_API mdac_status mdac_csda_set_imag(mdac_csda* v, const mdac_sda* x);
/* The operator<< text (the real part's table, then the imaginary part's). */
MDAC_API mdac_status mdac_csda_to_string(const mdac_csda* v, char* buf, size_t cap, size_t* n);

/* ---- CSDA: arithmetic ----------------------------------------------------------
 * Operands: c a CSDA, s an SDA, e an Expr, d a double, z a complex (re, im),
 * named as for CNDA: mdac_csda_<op> (c, c), _s/s<op>, _e/e<op>, _d/d<op>,
 * _z/z<op>; mdac_sda_<op>_z and mdac_sda_z<op> (s, z)/(z, s) give a CSDA too.
 * _into forms write into the existing CSDA out (same env; out may be an
 * operand). */

#define MDAC_CSDA_BINOP_DECL(op)                                                                        \
    MDAC_API mdac_status mdac_csda_##op(const mdac_csda* a, const mdac_csda* b, mdac_csda** out);       \
    MDAC_API mdac_status mdac_csda_##op##_s(const mdac_csda* a, const mdac_sda* b, mdac_csda** out);    \
    MDAC_API mdac_status mdac_csda_s##op(const mdac_sda* a, const mdac_csda* b, mdac_csda** out);       \
    MDAC_API mdac_status mdac_csda_##op##_e(const mdac_csda* a, const mdac_expr* x, mdac_csda** out);   \
    MDAC_API mdac_status mdac_csda_e##op(const mdac_expr* x, const mdac_csda* a, mdac_csda** out);      \
    MDAC_API mdac_status mdac_csda_##op##_d(const mdac_csda* a, double x, mdac_csda** out);             \
    MDAC_API mdac_status mdac_csda_d##op(double x, const mdac_csda* a, mdac_csda** out);                \
    MDAC_API mdac_status mdac_csda_##op##_z(const mdac_csda* a, double re, double im, mdac_csda** out); \
    MDAC_API mdac_status mdac_csda_z##op(double re, double im, const mdac_csda* a, mdac_csda** out);    \
    MDAC_API mdac_status mdac_sda_##op##_z(const mdac_sda* a, double re, double im, mdac_csda** out);   \
    MDAC_API mdac_status mdac_sda_z##op(double re, double im, const mdac_sda* a, mdac_csda** out);      \
    MDAC_API mdac_status mdac_csda_##op##_into(mdac_csda* out, const mdac_csda* a, const mdac_csda* b); \
    MDAC_API mdac_status mdac_csda_##op##_s_into(mdac_csda* out, const mdac_csda* a, const mdac_sda* b); \
    MDAC_API mdac_status mdac_csda_s##op##_into(mdac_csda* out, const mdac_sda* a, const mdac_csda* b); \
    MDAC_API mdac_status mdac_csda_##op##_e_into(mdac_csda* out, const mdac_csda* a, const mdac_expr* x); \
    MDAC_API mdac_status mdac_csda_e##op##_into(mdac_csda* out, const mdac_expr* x, const mdac_csda* a); \
    MDAC_API mdac_status mdac_csda_##op##_d_into(mdac_csda* out, const mdac_csda* a, double x);         \
    MDAC_API mdac_status mdac_csda_d##op##_into(mdac_csda* out, double x, const mdac_csda* a);          \
    MDAC_API mdac_status mdac_csda_##op##_z_into(mdac_csda* out, const mdac_csda* a, double re, double im); \
    MDAC_API mdac_status mdac_csda_z##op##_into(mdac_csda* out, double re, double im, const mdac_csda* a); \
    MDAC_API mdac_status mdac_sda_##op##_z_into(mdac_csda* out, const mdac_sda* a, double re, double im); \
    MDAC_API mdac_status mdac_sda_z##op##_into(mdac_csda* out, double re, double im, const mdac_sda* a);

MDAC_CSDA_BINOP_DECL(add)
MDAC_CSDA_BINOP_DECL(sub)
MDAC_CSDA_BINOP_DECL(mul)
MDAC_CSDA_BINOP_DECL(div)

MDAC_API mdac_status mdac_csda_neg(const mdac_csda* a, mdac_csda** out);
MDAC_API mdac_status mdac_csda_neg_into(mdac_csda* out, const mdac_csda* a);
MDAC_API mdac_status mdac_csda_pow_i(const mdac_csda* a, int n, mdac_csda** out);
MDAC_API mdac_status mdac_csda_pow_i_into(mdac_csda* out, const mdac_csda* a, int n);
MDAC_API mdac_status mdac_csda_pow_d(const mdac_csda* a, double x, mdac_csda** out);
MDAC_API mdac_status mdac_csda_pow_d_into(mdac_csda* out, const mdac_csda* a, double x);

/* ---- CSDA: math functions (the CNDA set) ---------------------------------------- */

#define MDAC_CSDA_FUNC_DECL(f)                                                   \
    MDAC_API mdac_status mdac_csda_##f(const mdac_csda* a, mdac_csda** out);     \
    MDAC_API mdac_status mdac_csda_##f##_into(mdac_csda* out, const mdac_csda* a);

MDAC_CSDA_FUNC_DECL(sqrt)
MDAC_CSDA_FUNC_DECL(exp)
MDAC_CSDA_FUNC_DECL(log)
MDAC_CSDA_FUNC_DECL(asin)
MDAC_CSDA_FUNC_DECL(acos)
MDAC_CSDA_FUNC_DECL(atan)
MDAC_CSDA_FUNC_DECL(asinh)
MDAC_CSDA_FUNC_DECL(acosh)
MDAC_CSDA_FUNC_DECL(atanh)

/* abs(CSDA): the modulus sqrt(re^2 + im^2), an SDA. */
MDAC_API mdac_status mdac_csda_abs(const mdac_csda* a, mdac_sda** out);

/* The CNDA of v with vals[j] for the symbol keys[j] (j < n), as mdac_sda_evaluate. */
MDAC_API mdac_status mdac_csda_evaluate(const mdac_csda* v, size_t n, const mdac_expr* const* keys,
                                        const double* vals, mdac_cnda** out);

/* ---- CSDA lists and composition (as CNDA) ---------------------------------------- */

MDAC_API mdac_status mdac_csdalist_new(mdac_csdalist** out);
MDAC_API mdac_status mdac_csdalist_from(const mdac_csda* const* vs, size_t n, mdac_csdalist** out);
MDAC_API void mdac_csdalist_free(mdac_csdalist* l);
MDAC_API size_t mdac_csdalist_length(const mdac_csdalist* l);
MDAC_API mdac_env* mdac_csdalist_env(const mdac_csdalist* l);
MDAC_API mdac_status mdac_csdalist_get(const mdac_csdalist* l, size_t i, mdac_csda** out);
MDAC_API mdac_status mdac_csdalist_set(mdac_csdalist* l, size_t i, const mdac_csda* v);
MDAC_API mdac_status mdac_csdalist_push(mdac_csdalist* l, const mdac_csda* v);

/* An SDA map at CSDA arguments. */
MDAC_API mdac_status mdac_sdalist_compose_c(const mdac_sdalist* m, const mdac_csdalist* v,
                                            mdac_csdalist** out);
/* A CSDA map at CSDA arguments. */
MDAC_API mdac_status mdac_csdalist_compose(const mdac_csdalist* m, const mdac_csdalist* v,
                                           mdac_csdalist** out);
/* A CSDA map at SDA arguments. */
MDAC_API mdac_status mdac_csdalist_compose_s(const mdac_csdalist* m, const mdac_sdalist* v,
                                             mdac_csdalist** out);

#ifdef __cplusplus
}
#endif

#endif /* MIRADAC_H */
