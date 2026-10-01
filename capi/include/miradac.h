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

#ifdef __cplusplus
}
#endif

#endif /* MIRADAC_H */
