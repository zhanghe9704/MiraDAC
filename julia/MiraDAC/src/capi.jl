# capi.jl — the C API (capi/include/miradac.h) through ccall, status checking, the deferred
# free queues and the allocation retry (DEVELOPMENT_PLAN_JULIA.md A.3, A.5).

const Handle = Ptr{Cvoid}

# ---- Status codes and exceptions ------------------------------------------------------------

const MDAC_OK = Cint(0)
const MDAC_ERR_ENV = Cint(1)
const MDAC_ERR_VALUE = Cint(2)
const MDAC_ERR_INDEX = Cint(3)
const MDAC_ERR_POOL = Cint(4)
const MDAC_ERR_RUNTIME = Cint(5)
const MDAC_ERR_UNSUPPORTED = Cint(6)

"""Supertype of the exceptions MiraDAC throws besides the Base ones."""
abstract type MiraDACError <: Exception end

"""A DA environment mismatch, a cleared environment, or no environment."""
struct EnvError <: MiraDACError
    msg::String
end

"""The environment's pool has no free vector, even after a garbage collection."""
struct PoolExhaustedError <: MiraDACError
    msg::String
end

Base.showerror(io::IO, e::MiraDACError) = print(io, nameof(typeof(e)), ": ", e.msg)

mdac_last_error() = unsafe_string(ccall((:mdac_last_error, libmiradac), Cstring, ()))

"""
    check(st)
    check(st, x)

Throw the exception of a failed C API call (A.6). With `x`, `MDAC_ERR_VALUE` becomes a
`DomainError` of `x` (math functions); otherwise an `ArgumentError`.
"""
@inline check(st::Cint) = st == MDAC_OK ? nothing : throw_status(st, nothing)
@inline check(st::Cint, x) = st == MDAC_OK ? nothing : throw_status(st, x)

@noinline function throw_status(st::Cint, x)
    msg = mdac_last_error()
    st == MDAC_ERR_ENV && throw(EnvError(msg))
    st == MDAC_ERR_VALUE && throw(x === nothing ? ArgumentError(msg) : DomainError(x, msg))
    st == MDAC_ERR_INDEX && throw(BoundsError())
    st == MDAC_ERR_POOL && throw(PoolExhaustedError(msg))
    throw(ErrorException(msg))
end

# ---- Versions and envs ----------------------------------------------------------------------

mdac_abi_version() = ccall((:mdac_abi_version, libmiradac), Cint, ())
mdac_version() = unsafe_string(ccall((:mdac_version, libmiradac), Cstring, ()))

mdac_init(order, nvars, poolsize, table) =
    ccall((:mdac_init, libmiradac), Cint, (Cuint, Cuint, Cuint, Cint), order, nvars, poolsize, table)
mdac_clear() = ccall((:mdac_clear, libmiradac), Cint, ())

mdac_env_current(out) = ccall((:mdac_env_current, libmiradac), Cint, (Ptr{Handle},), out)
mdac_env_default(out) = ccall((:mdac_env_default, libmiradac), Cint, (Ptr{Handle},), out)
mdac_env_make(order, nvars, poolsize, table, out) =
    ccall((:mdac_env_make, libmiradac), Cint, (Cuint, Cuint, Cuint, Cint, Ptr{Handle}),
          order, nvars, poolsize, table, out)
mdac_env_select(e) = ccall((:mdac_env_select, libmiradac), Cint, (Handle,), e)
mdac_env_exchange(e, prev) =
    ccall((:mdac_env_exchange, libmiradac), Cint, (Handle, Ptr{Handle}), e, prev)
mdac_env_close(e) = ccall((:mdac_env_close, libmiradac), Cint, (Handle,), e)

for q in (:order, :max_order, :nvars, :full_length, :poolsize, :count, :remain)
    f = Symbol(:mdac_env_, q)
    @eval $f(e, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Cuint}), e, out)
end
mdac_env_retired(e) = ccall((:mdac_env_retired, libmiradac), Cint, (Handle,), e)
mdac_env_change_order(e, n, ok) =
    ccall((:mdac_env_change_order, libmiradac), Cint, (Handle, Cuint, Ptr{Cint}), e, n, ok)
mdac_env_restore_order(e) = ccall((:mdac_env_restore_order, libmiradac), Cint, (Handle,), e)

mdac_get_eps() = ccall((:mdac_get_eps, libmiradac), Cdouble, ())
mdac_set_eps(x) = ccall((:mdac_set_eps, libmiradac), Cint, (Cdouble,), x)

# ---- NDA: lifecycle and inspection ----------------------------------------------------------

mdac_nda_new(e, x, out) =
    ccall((:mdac_nda_new, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), e, x, out)
mdac_nda_from_coeffs(e, c, n, out) =
    ccall((:mdac_nda_from_coeffs, libmiradac), Cint, (Handle, Ptr{Cdouble}, Csize_t, Ptr{Handle}),
          e, c, n, out)
mdac_nda_var(e, i, out) =
    ccall((:mdac_nda_var, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), e, i, out)
mdac_nda_copy(v, out) = ccall((:mdac_nda_copy, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_nda_free(v) = ccall((:mdac_nda_free, libmiradac), Cvoid, (Handle,), v)
mdac_nda_env(v) = ccall((:mdac_nda_env, libmiradac), Handle, (Handle,), v)

mdac_nda_con(v, out) = ccall((:mdac_nda_con, libmiradac), Cint, (Handle, Ptr{Cdouble}), v, out)
mdac_nda_set_con(v, x) = ccall((:mdac_nda_set_con, libmiradac), Cint, (Handle, Cdouble), v, x)
mdac_nda_length(v, out) = ccall((:mdac_nda_length, libmiradac), Cint, (Handle, Ptr{Csize_t}), v, out)
mdac_nda_nterms(v, out) = ccall((:mdac_nda_nterms, libmiradac), Cint, (Handle, Ptr{Csize_t}), v, out)
mdac_nda_norm(v, out) = ccall((:mdac_nda_norm, libmiradac), Cint, (Handle, Ptr{Cdouble}), v, out)
mdac_nda_coeffs(v, buf, cap, n) =
    ccall((:mdac_nda_coeffs, libmiradac), Cint, (Handle, Ptr{Cdouble}, Csize_t, Ptr{Csize_t}),
          v, buf, cap, n)
mdac_nda_coeff(v, exps, k, out) =
    ccall((:mdac_nda_coeff, libmiradac), Cint, (Handle, Ptr{Cint}, Csize_t, Ptr{Cdouble}),
          v, exps, k, out)
mdac_nda_set_coeff(v, exps, k, x) =
    ccall((:mdac_nda_set_coeff, libmiradac), Cint, (Handle, Ptr{Cint}, Csize_t, Cdouble),
          v, exps, k, x)
mdac_nda_index_term(v, i, exps, out) =
    ccall((:mdac_nda_index_term, libmiradac), Cint, (Handle, Csize_t, Ptr{Cint}, Ptr{Cdouble}),
          v, i, exps, out)
mdac_nda_iszero(v, eps, out) =
    ccall((:mdac_nda_iszero, libmiradac), Cint, (Handle, Cdouble, Ptr{Cint}), v, eps, out)
mdac_nda_clean(v, eps) = ccall((:mdac_nda_clean, libmiradac), Cint, (Handle, Cdouble), v, eps)
mdac_nda_reset(v) = ccall((:mdac_nda_reset, libmiradac), Cint, (Handle,), v)
mdac_nda_to_string(v, buf, cap, n) =
    ccall((:mdac_nda_to_string, libmiradac), Cint, (Handle, Ptr{UInt8}, Csize_t, Ptr{Csize_t}),
          v, buf, cap, n)

# ---- NDA: arithmetic ------------------------------------------------------------------------

for op in (:add, :sub, :mul, :div)
    f, fd, df = Symbol(:mdac_nda_, op), Symbol(:mdac_nda_, op, :_d), Symbol(:mdac_nda_d, op)
    fi, fdi, dfi = Symbol(f, :_into), Symbol(fd, :_into), Symbol(df, :_into)
    @eval begin
        $f(a, b, out) = ccall(($(QuoteNode(f)), libmiradac), Cint,
                              (Handle, Handle, Ptr{Handle}), a, b, out)
        $fd(a, x, out) = ccall(($(QuoteNode(fd)), libmiradac), Cint,
                               (Handle, Cdouble, Ptr{Handle}), a, x, out)
        $df(x, a, out) = ccall(($(QuoteNode(df)), libmiradac), Cint,
                               (Cdouble, Handle, Ptr{Handle}), x, a, out)
        $fi(out, a, b) = ccall(($(QuoteNode(fi)), libmiradac), Cint,
                               (Handle, Handle, Handle), out, a, b)
        $fdi(out, a, x) = ccall(($(QuoteNode(fdi)), libmiradac), Cint,
                                (Handle, Handle, Cdouble), out, a, x)
        $dfi(out, x, a) = ccall(($(QuoteNode(dfi)), libmiradac), Cint,
                                (Handle, Cdouble, Handle), out, x, a)
    end
end

mdac_nda_neg(a, out) = ccall((:mdac_nda_neg, libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
mdac_nda_neg_into(out, a) = ccall((:mdac_nda_neg_into, libmiradac), Cint, (Handle, Handle), out, a)
mdac_nda_pow_i(a, n, out) =
    ccall((:mdac_nda_pow_i, libmiradac), Cint, (Handle, Cint, Ptr{Handle}), a, n, out)
mdac_nda_pow_i_into(out, a, n) =
    ccall((:mdac_nda_pow_i_into, libmiradac), Cint, (Handle, Handle, Cint), out, a, n)
mdac_nda_pow_d(a, x, out) =
    ccall((:mdac_nda_pow_d, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), a, x, out)
mdac_nda_pow_d_into(out, a, x) =
    ccall((:mdac_nda_pow_d_into, libmiradac), Cint, (Handle, Handle, Cdouble), out, a, x)

# ---- NDA: math functions --------------------------------------------------------------------

const NDA_FUNCS = (:sqrt, :exp, :log, :sin, :cos, :tan, :asin, :acos, :atan,
                   :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :erf)

for fn in NDA_FUNCS
    f, fi = Symbol(:mdac_nda_, fn), Symbol(:mdac_nda_, fn, :_into)
    @eval begin
        $f(a, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
        $fi(out, a) = ccall(($(QuoteNode(fi)), libmiradac), Cint, (Handle, Handle), out, a)
    end
end

mdac_nda_abs(a, out) = ccall((:mdac_nda_abs, libmiradac), Cint, (Handle, Ptr{Cdouble}), a, out)

# ---- NDA lists and algorithms ---------------------------------------------------------------

mdac_ndalist_new(out) = ccall((:mdac_ndalist_new, libmiradac), Cint, (Ptr{Handle},), out)
mdac_ndalist_from(vs, n, out) =
    ccall((:mdac_ndalist_from, libmiradac), Cint, (Ptr{Handle}, Csize_t, Ptr{Handle}), vs, n, out)
mdac_ndalist_free(l) = ccall((:mdac_ndalist_free, libmiradac), Cvoid, (Handle,), l)
mdac_ndalist_length(l) = ccall((:mdac_ndalist_length, libmiradac), Csize_t, (Handle,), l)
mdac_ndalist_env(l) = ccall((:mdac_ndalist_env, libmiradac), Handle, (Handle,), l)
mdac_ndalist_get(l, i, out) =
    ccall((:mdac_ndalist_get, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}), l, i, out)
mdac_ndalist_set(l, i, v) =
    ccall((:mdac_ndalist_set, libmiradac), Cint, (Handle, Csize_t, Handle), l, i, v)
mdac_ndalist_push(l, v) = ccall((:mdac_ndalist_push, libmiradac), Cint, (Handle, Handle), l, v)

mdac_nda_der(v, i, out) =
    ccall((:mdac_nda_der, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), v, i, out)
mdac_nda_integ(v, i, out) =
    ccall((:mdac_nda_integ, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), v, i, out)
mdac_nda_substitute_d(v, i, x, out) =
    ccall((:mdac_nda_substitute_d, libmiradac), Cint, (Handle, Cuint, Cdouble, Ptr{Handle}),
          v, i, x, out)
mdac_nda_substitute(v, i, x, out) =
    ccall((:mdac_nda_substitute, libmiradac), Cint, (Handle, Cuint, Handle, Ptr{Handle}),
          v, i, x, out)
mdac_nda_substitute_multi(v, ids, k, xs, out) =
    ccall((:mdac_nda_substitute_multi, libmiradac), Cint,
          (Handle, Ptr{Cuint}, Csize_t, Handle, Ptr{Handle}), v, ids, k, xs, out)
mdac_ndalist_substitute(m, ids, k, xs, out) =
    ccall((:mdac_ndalist_substitute, libmiradac), Cint,
          (Handle, Ptr{Cuint}, Csize_t, Handle, Ptr{Handle}), m, ids, k, xs, out)
mdac_ndalist_compose(m, v, out) =
    ccall((:mdac_ndalist_compose, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), m, v, out)
mdac_ndalist_compose_d(m, pt, n, out) =
    ccall((:mdac_ndalist_compose_d, libmiradac), Cint,
          (Handle, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}), m, pt, n, out)
mdac_ndalist_compose_z(m, pt, n, out) =
    ccall((:mdac_ndalist_compose_z, libmiradac), Cint,
          (Handle, Ptr{ComplexF64}, Csize_t, Ptr{ComplexF64}), m, pt, n, out)
mdac_ndalist_inv_map(m, dim, out) =
    ccall((:mdac_ndalist_inv_map, libmiradac), Cint, (Handle, Cint, Ptr{Handle}), m, dim, out)
mdac_ndalist_evaluate_map(m, pts, npts, out) =
    ccall((:mdac_ndalist_evaluate_map, libmiradac), Cint,
          (Handle, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}), m, pts, npts, out)
mdac_nda_eval(v, pt, n, out) =
    ccall((:mdac_nda_eval, libmiradac), Cint, (Handle, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}),
          v, pt, n, out)
mdac_exponents(e, buf, cap, n) =
    ccall((:mdac_exponents, libmiradac), Cint, (Handle, Ptr{Cint}, Csize_t, Ptr{Csize_t}),
          e, buf, cap, n)

# ---- CNDA -----------------------------------------------------------------------------------

mdac_cnda_new(re, im, out) =
    ccall((:mdac_cnda_new, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), re, im, out)
mdac_cnda_new_z(e, re, im, out) =
    ccall((:mdac_cnda_new_z, libmiradac), Cint, (Handle, Cdouble, Cdouble, Ptr{Handle}),
          e, re, im, out)
mdac_cnda_copy(v, out) = ccall((:mdac_cnda_copy, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_cnda_free(v) = ccall((:mdac_cnda_free, libmiradac), Cvoid, (Handle,), v)
mdac_cnda_env(v) = ccall((:mdac_cnda_env, libmiradac), Handle, (Handle,), v)
mdac_cnda_real(v, out) = ccall((:mdac_cnda_real, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_cnda_imag(v, out) = ccall((:mdac_cnda_imag, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_cnda_set_real(v, x) = ccall((:mdac_cnda_set_real, libmiradac), Cint, (Handle, Handle), v, x)
mdac_cnda_set_imag(v, x) = ccall((:mdac_cnda_set_imag, libmiradac), Cint, (Handle, Handle), v, x)
mdac_cnda_to_string(v, buf, cap, n) =
    ccall((:mdac_cnda_to_string, libmiradac), Cint, (Handle, Ptr{UInt8}, Csize_t, Ptr{Csize_t}),
          v, buf, cap, n)

# The operand shapes of a CNDA operator (miradac.h): C function name => argument types; a
# complex operand is two doubles.
cnda_op_shapes(op) = (
    Symbol(:mdac_cnda_, op) => (Handle, Handle),
    Symbol(:mdac_cnda_, op, :_n) => (Handle, Handle),
    Symbol(:mdac_cnda_n, op) => (Handle, Handle),
    Symbol(:mdac_cnda_, op, :_d) => (Handle, Cdouble),
    Symbol(:mdac_cnda_d, op) => (Cdouble, Handle),
    Symbol(:mdac_cnda_, op, :_z) => (Handle, Cdouble, Cdouble),
    Symbol(:mdac_cnda_z, op) => (Cdouble, Cdouble, Handle),
    Symbol(:mdac_nda_, op, :_z) => (Handle, Cdouble, Cdouble),
    Symbol(:mdac_nda_z, op) => (Cdouble, Cdouble, Handle))

for op in (:add, :sub, :mul, :div), (f, types) in cnda_op_shapes(op)
    fi = Symbol(f, :_into)
    args = [Symbol(:x, k) for k in eachindex(types)]
    @eval begin
        $f($(args...), out) = ccall(($(QuoteNode(f)), libmiradac), Cint,
                                    ($(types...), Ptr{Handle}), $(args...), out)
        $fi(out, $(args...)) = ccall(($(QuoteNode(fi)), libmiradac), Cint,
                                     (Handle, $(types...)), out, $(args...))
    end
end

mdac_cnda_neg(a, out) = ccall((:mdac_cnda_neg, libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
mdac_cnda_neg_into(out, a) = ccall((:mdac_cnda_neg_into, libmiradac), Cint, (Handle, Handle), out, a)
mdac_cnda_pow_i(a, n, out) =
    ccall((:mdac_cnda_pow_i, libmiradac), Cint, (Handle, Cint, Ptr{Handle}), a, n, out)
mdac_cnda_pow_i_into(out, a, n) =
    ccall((:mdac_cnda_pow_i_into, libmiradac), Cint, (Handle, Handle, Cint), out, a, n)
mdac_cnda_pow_d(a, x, out) =
    ccall((:mdac_cnda_pow_d, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), a, x, out)
mdac_cnda_pow_d_into(out, a, x) =
    ccall((:mdac_cnda_pow_d_into, libmiradac), Cint, (Handle, Handle, Cdouble), out, a, x)

# C++ has no complex sin, cos, tan, sinh, cosh, tanh or erf.
const CNDA_FUNCS = (:sqrt, :exp, :log, :asin, :acos, :atan, :asinh, :acosh, :atanh)

for fn in CNDA_FUNCS
    f, fi = Symbol(:mdac_cnda_, fn), Symbol(:mdac_cnda_, fn, :_into)
    @eval begin
        $f(a, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
        $fi(out, a) = ccall(($(QuoteNode(fi)), libmiradac), Cint, (Handle, Handle), out, a)
    end
end

mdac_cnda_abs(a, out) = ccall((:mdac_cnda_abs, libmiradac), Cint, (Handle, Ptr{Cdouble}), a, out)

mdac_cndalist_new(out) = ccall((:mdac_cndalist_new, libmiradac), Cint, (Ptr{Handle},), out)
mdac_cndalist_from(vs, n, out) =
    ccall((:mdac_cndalist_from, libmiradac), Cint, (Ptr{Handle}, Csize_t, Ptr{Handle}), vs, n, out)
mdac_cndalist_free(l) = ccall((:mdac_cndalist_free, libmiradac), Cvoid, (Handle,), l)
mdac_cndalist_length(l) = ccall((:mdac_cndalist_length, libmiradac), Csize_t, (Handle,), l)
mdac_cndalist_env(l) = ccall((:mdac_cndalist_env, libmiradac), Handle, (Handle,), l)
mdac_cndalist_get(l, i, out) =
    ccall((:mdac_cndalist_get, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}), l, i, out)
mdac_cndalist_set(l, i, v) =
    ccall((:mdac_cndalist_set, libmiradac), Cint, (Handle, Csize_t, Handle), l, i, v)
mdac_cndalist_push(l, v) = ccall((:mdac_cndalist_push, libmiradac), Cint, (Handle, Handle), l, v)

for f in (:mdac_ndalist_compose_c, :mdac_cndalist_compose, :mdac_cndalist_compose_n)
    @eval $f(m, v, out) = ccall(($(QuoteNode(f)), libmiradac), Cint,
                                (Handle, Handle, Ptr{Handle}), m, v, out)
end

# ---- Symbolic: Expr --------------------------------------------------------------------------

mdac_has_symbolic() = ccall((:mdac_has_symbolic, libmiradac), Cint, ())

mdac_expr_new_d(x, out) = ccall((:mdac_expr_new_d, libmiradac), Cint, (Cdouble, Ptr{Handle}), x, out)
mdac_expr_new_i(x, out) = ccall((:mdac_expr_new_i, libmiradac), Cint, (Int64, Ptr{Handle}), x, out)
mdac_expr_parse(s, out) = ccall((:mdac_expr_parse, libmiradac), Cint, (Cstring, Ptr{Handle}), s, out)
mdac_expr_symbol(s, out) = ccall((:mdac_expr_symbol, libmiradac), Cint, (Cstring, Ptr{Handle}), s, out)
mdac_expr_copy(x, out) = ccall((:mdac_expr_copy, libmiradac), Cint, (Handle, Ptr{Handle}), x, out)
mdac_expr_free(x) = ccall((:mdac_expr_free, libmiradac), Cvoid, (Handle,), x)
mdac_expr_to_string(x, buf, cap, n) =
    ccall((:mdac_expr_to_string, libmiradac), Cint, (Handle, Ptr{UInt8}, Csize_t, Ptr{Csize_t}),
          x, buf, cap, n)
mdac_expr_to_double(x, out) =
    ccall((:mdac_expr_to_double, libmiradac), Cint, (Handle, Ptr{Cdouble}), x, out)

for op in (:add, :sub, :mul, :div, :pow)
    f, fd, df = Symbol(:mdac_expr_, op), Symbol(:mdac_expr_, op, :_d), Symbol(:mdac_expr_d, op)
    @eval begin
        $f(a, b, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Handle, Ptr{Handle}), a, b, out)
        $fd(a, x, out) = ccall(($(QuoteNode(fd)), libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), a, x, out)
        $df(x, a, out) = ccall(($(QuoteNode(df)), libmiradac), Cint, (Cdouble, Handle, Ptr{Handle}), x, a, out)
    end
end

mdac_expr_neg(a, out) = ccall((:mdac_expr_neg, libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
mdac_expr_eq(a, b, out) = ccall((:mdac_expr_eq, libmiradac), Cint, (Handle, Handle, Ptr{Cint}), a, b, out)
mdac_expr_hash(a, out) = ccall((:mdac_expr_hash, libmiradac), Cint, (Handle, Ptr{UInt64}), a, out)
mdac_expr_is_zero(a, out) = ccall((:mdac_expr_is_zero, libmiradac), Cint, (Handle, Ptr{Cint}), a, out)
mdac_expr_subs(x, n, keys, vals, out) =
    ccall((:mdac_expr_subs, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}, Ptr{Handle}, Ptr{Handle}),
          x, n, keys, vals, out)
for f in (:mdac_expr_expand, :mdac_expr_simplify)
    @eval $f(x, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), x, out)
end
mdac_expr_diff(x, s, out) =
    ccall((:mdac_expr_diff, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), x, s, out)
mdac_expr_free_symbols(x, buf, cap, n) =
    ccall((:mdac_expr_free_symbols, libmiradac), Cint, (Handle, Ptr{Handle}, Csize_t, Ptr{Csize_t}),
          x, buf, cap, n)

# ---- Symbolic: SDA ---------------------------------------------------------------------------

mdac_sda_new(e, x, out) = ccall((:mdac_sda_new, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), e, x, out)
mdac_sda_new_d(e, x, out) =
    ccall((:mdac_sda_new_d, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), e, x, out)
mdac_sda_var(e, i, out) = ccall((:mdac_sda_var, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), e, i, out)
mdac_sda_promote(v, out) = ccall((:mdac_sda_promote, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_sda_copy(v, out) = ccall((:mdac_sda_copy, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_sda_free(v) = ccall((:mdac_sda_free, libmiradac), Cvoid, (Handle,), v)
mdac_sda_env(v) = ccall((:mdac_sda_env, libmiradac), Handle, (Handle,), v)
mdac_sda_con(v, out) = ccall((:mdac_sda_con, libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
mdac_sda_length(v, out) = ccall((:mdac_sda_length, libmiradac), Cint, (Handle, Ptr{Csize_t}), v, out)
mdac_sda_nterms(v, out) = ccall((:mdac_sda_nterms, libmiradac), Cint, (Handle, Ptr{Csize_t}), v, out)
mdac_sda_coeffs(v, buf, cap, n) =
    ccall((:mdac_sda_coeffs, libmiradac), Cint, (Handle, Ptr{Handle}, Csize_t, Ptr{Csize_t}), v, buf, cap, n)
mdac_sda_coeff(v, exps, k, out) =
    ccall((:mdac_sda_coeff, libmiradac), Cint, (Handle, Ptr{Cint}, Csize_t, Ptr{Handle}), v, exps, k, out)
mdac_sda_index_term(v, i, exps, out) =
    ccall((:mdac_sda_index_term, libmiradac), Cint, (Handle, Csize_t, Ptr{Cint}, Ptr{Handle}),
          v, i, exps, out)
mdac_sda_iszero(v, out) = ccall((:mdac_sda_iszero, libmiradac), Cint, (Handle, Ptr{Cint}), v, out)
mdac_sda_to_string(v, buf, cap, n) =
    ccall((:mdac_sda_to_string, libmiradac), Cint, (Handle, Ptr{UInt8}, Csize_t, Ptr{Csize_t}),
          v, buf, cap, n)

# The operand shapes of an SDA operator (miradac.h), as cnda_op_shapes.
sda_op_shapes(op) = (
    Symbol(:mdac_sda_, op) => (Handle, Handle),
    Symbol(:mdac_sda_, op, :_n) => (Handle, Handle),
    Symbol(:mdac_sda_n, op) => (Handle, Handle),
    Symbol(:mdac_sda_, op, :_e) => (Handle, Handle),
    Symbol(:mdac_sda_e, op) => (Handle, Handle),
    Symbol(:mdac_sda_, op, :_d) => (Handle, Cdouble),
    Symbol(:mdac_sda_d, op) => (Cdouble, Handle),
    Symbol(:mdac_nda_, op, :_e) => (Handle, Handle),
    Symbol(:mdac_nda_e, op) => (Handle, Handle))

for op in (:add, :sub, :mul, :div), (f, types) in sda_op_shapes(op)
    fi = Symbol(f, :_into)
    @eval begin
        $f(x1, x2, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, ($(types...), Ptr{Handle}), x1, x2, out)
        $fi(out, x1, x2) = ccall(($(QuoteNode(fi)), libmiradac), Cint, (Handle, $(types...)), out, x1, x2)
    end
end

mdac_sda_neg(a, out) = ccall((:mdac_sda_neg, libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
mdac_sda_neg_into(out, a) = ccall((:mdac_sda_neg_into, libmiradac), Cint, (Handle, Handle), out, a)
mdac_sda_pow_i(a, n, out) =
    ccall((:mdac_sda_pow_i, libmiradac), Cint, (Handle, Cint, Ptr{Handle}), a, n, out)
mdac_sda_pow_i_into(out, a, n) =
    ccall((:mdac_sda_pow_i_into, libmiradac), Cint, (Handle, Handle, Cint), out, a, n)
mdac_sda_pow_d(a, x, out) =
    ccall((:mdac_sda_pow_d, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), a, x, out)
mdac_sda_pow_d_into(out, a, x) =
    ccall((:mdac_sda_pow_d_into, libmiradac), Cint, (Handle, Handle, Cdouble), out, a, x)

# C++ has no SDA asinh, acosh or atanh.
const SDA_FUNCS = (:sqrt, :exp, :log, :sin, :cos, :tan, :asin, :acos, :atan,
                   :sinh, :cosh, :tanh, :erf)

for fn in SDA_FUNCS
    f, fi = Symbol(:mdac_sda_, fn), Symbol(:mdac_sda_, fn, :_into)
    @eval begin
        $f(a, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
        $fi(out, a) = ccall(($(QuoteNode(fi)), libmiradac), Cint, (Handle, Handle), out, a)
    end
end

for f in (:mdac_sda_simplify, :mdac_sda_expand)
    @eval $f(v, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
end
mdac_sda_subs(v, n, keys, vals, out) =
    ccall((:mdac_sda_subs, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}, Ptr{Handle}, Ptr{Handle}),
          v, n, keys, vals, out)
mdac_sda_evaluate(v, n, keys, vals, out) =
    ccall((:mdac_sda_evaluate, libmiradac), Cint,
          (Handle, Csize_t, Ptr{Handle}, Ptr{Cdouble}, Ptr{Handle}), v, n, keys, vals, out)

mdac_sdalist_new(out) = ccall((:mdac_sdalist_new, libmiradac), Cint, (Ptr{Handle},), out)
mdac_sdalist_from(vs, n, out) =
    ccall((:mdac_sdalist_from, libmiradac), Cint, (Ptr{Handle}, Csize_t, Ptr{Handle}), vs, n, out)
mdac_sdalist_free(l) = ccall((:mdac_sdalist_free, libmiradac), Cvoid, (Handle,), l)
mdac_sdalist_length(l) = ccall((:mdac_sdalist_length, libmiradac), Csize_t, (Handle,), l)
mdac_sdalist_env(l) = ccall((:mdac_sdalist_env, libmiradac), Handle, (Handle,), l)
mdac_sdalist_get(l, i, out) =
    ccall((:mdac_sdalist_get, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}), l, i, out)
mdac_sdalist_set(l, i, v) =
    ccall((:mdac_sdalist_set, libmiradac), Cint, (Handle, Csize_t, Handle), l, i, v)
mdac_sdalist_push(l, v) = ccall((:mdac_sdalist_push, libmiradac), Cint, (Handle, Handle), l, v)

mdac_sda_der(v, i, out) = ccall((:mdac_sda_der, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), v, i, out)
mdac_sda_integ(v, i, out) =
    ccall((:mdac_sda_integ, libmiradac), Cint, (Handle, Cuint, Ptr{Handle}), v, i, out)
mdac_sda_substitute_d(v, i, x, out) =
    ccall((:mdac_sda_substitute_d, libmiradac), Cint, (Handle, Cuint, Cdouble, Ptr{Handle}), v, i, x, out)
mdac_sda_substitute(v, i, x, out) =
    ccall((:mdac_sda_substitute, libmiradac), Cint, (Handle, Cuint, Handle, Ptr{Handle}), v, i, x, out)
mdac_sda_substitute_multi(v, ids, k, xs, out) =
    ccall((:mdac_sda_substitute_multi, libmiradac), Cint,
          (Handle, Ptr{Cuint}, Csize_t, Handle, Ptr{Handle}), v, ids, k, xs, out)
mdac_sdalist_substitute(m, ids, k, xs, out) =
    ccall((:mdac_sdalist_substitute, libmiradac), Cint,
          (Handle, Ptr{Cuint}, Csize_t, Handle, Ptr{Handle}), m, ids, k, xs, out)
mdac_sdalist_compose(m, v, out) =
    ccall((:mdac_sdalist_compose, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), m, v, out)
mdac_sdalist_compose_d(m, pt, n, out) =
    ccall((:mdac_sdalist_compose_d, libmiradac), Cint, (Handle, Ptr{Cdouble}, Csize_t, Ptr{Handle}),
          m, pt, n, out)

# ---- Symbolic: CSDA --------------------------------------------------------------------------

mdac_csda_new(re, im, out) =
    ccall((:mdac_csda_new, libmiradac), Cint, (Handle, Handle, Ptr{Handle}), re, im, out)
mdac_csda_new_z(e, re, im, out) =
    ccall((:mdac_csda_new_z, libmiradac), Cint, (Handle, Cdouble, Cdouble, Ptr{Handle}), e, re, im, out)
mdac_csda_free(v) = ccall((:mdac_csda_free, libmiradac), Cvoid, (Handle,), v)
mdac_csda_env(v) = ccall((:mdac_csda_env, libmiradac), Handle, (Handle,), v)
for f in (:mdac_csda_promote, :mdac_csda_copy, :mdac_csda_real, :mdac_csda_imag, :mdac_csda_neg,
          :mdac_csda_abs)
    @eval $f(v, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), v, out)
end
mdac_csda_set_real(v, x) = ccall((:mdac_csda_set_real, libmiradac), Cint, (Handle, Handle), v, x)
mdac_csda_set_imag(v, x) = ccall((:mdac_csda_set_imag, libmiradac), Cint, (Handle, Handle), v, x)
mdac_csda_to_string(v, buf, cap, n) =
    ccall((:mdac_csda_to_string, libmiradac), Cint, (Handle, Ptr{UInt8}, Csize_t, Ptr{Csize_t}),
          v, buf, cap, n)

# The operand shapes of a CSDA operator (miradac.h), as cnda_op_shapes.
csda_op_shapes(op) = (
    Symbol(:mdac_csda_, op) => (Handle, Handle),
    Symbol(:mdac_csda_, op, :_s) => (Handle, Handle),
    Symbol(:mdac_csda_s, op) => (Handle, Handle),
    Symbol(:mdac_csda_, op, :_e) => (Handle, Handle),
    Symbol(:mdac_csda_e, op) => (Handle, Handle),
    Symbol(:mdac_csda_, op, :_d) => (Handle, Cdouble),
    Symbol(:mdac_csda_d, op) => (Cdouble, Handle),
    Symbol(:mdac_csda_, op, :_z) => (Handle, Cdouble, Cdouble),
    Symbol(:mdac_csda_z, op) => (Cdouble, Cdouble, Handle),
    Symbol(:mdac_sda_, op, :_z) => (Handle, Cdouble, Cdouble),
    Symbol(:mdac_sda_z, op) => (Cdouble, Cdouble, Handle))

for op in (:add, :sub, :mul, :div), (f, types) in csda_op_shapes(op)
    fi = Symbol(f, :_into)
    args = [Symbol(:x, k) for k in eachindex(types)]
    @eval begin
        $f($(args...), out) = ccall(($(QuoteNode(f)), libmiradac), Cint,
                                    ($(types...), Ptr{Handle}), $(args...), out)
        $fi(out, $(args...)) = ccall(($(QuoteNode(fi)), libmiradac), Cint,
                                     (Handle, $(types...)), out, $(args...))
    end
end

mdac_csda_neg_into(out, a) = ccall((:mdac_csda_neg_into, libmiradac), Cint, (Handle, Handle), out, a)
mdac_csda_pow_i(a, n, out) =
    ccall((:mdac_csda_pow_i, libmiradac), Cint, (Handle, Cint, Ptr{Handle}), a, n, out)
mdac_csda_pow_i_into(out, a, n) =
    ccall((:mdac_csda_pow_i_into, libmiradac), Cint, (Handle, Handle, Cint), out, a, n)
mdac_csda_pow_d(a, x, out) =
    ccall((:mdac_csda_pow_d, libmiradac), Cint, (Handle, Cdouble, Ptr{Handle}), a, x, out)
mdac_csda_pow_d_into(out, a, x) =
    ccall((:mdac_csda_pow_d_into, libmiradac), Cint, (Handle, Handle, Cdouble), out, a, x)

for fn in CNDA_FUNCS
    f, fi = Symbol(:mdac_csda_, fn), Symbol(:mdac_csda_, fn, :_into)
    @eval begin
        $f(a, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Ptr{Handle}), a, out)
        $fi(out, a) = ccall(($(QuoteNode(fi)), libmiradac), Cint, (Handle, Handle), out, a)
    end
end

mdac_csda_evaluate(v, n, keys, vals, out) =
    ccall((:mdac_csda_evaluate, libmiradac), Cint,
          (Handle, Csize_t, Ptr{Handle}, Ptr{Cdouble}, Ptr{Handle}), v, n, keys, vals, out)

mdac_csdalist_new(out) = ccall((:mdac_csdalist_new, libmiradac), Cint, (Ptr{Handle},), out)
mdac_csdalist_from(vs, n, out) =
    ccall((:mdac_csdalist_from, libmiradac), Cint, (Ptr{Handle}, Csize_t, Ptr{Handle}), vs, n, out)
mdac_csdalist_free(l) = ccall((:mdac_csdalist_free, libmiradac), Cvoid, (Handle,), l)
mdac_csdalist_length(l) = ccall((:mdac_csdalist_length, libmiradac), Csize_t, (Handle,), l)
mdac_csdalist_env(l) = ccall((:mdac_csdalist_env, libmiradac), Handle, (Handle,), l)
mdac_csdalist_get(l, i, out) =
    ccall((:mdac_csdalist_get, libmiradac), Cint, (Handle, Csize_t, Ptr{Handle}), l, i, out)
mdac_csdalist_set(l, i, v) =
    ccall((:mdac_csdalist_set, libmiradac), Cint, (Handle, Csize_t, Handle), l, i, v)
mdac_csdalist_push(l, v) = ccall((:mdac_csdalist_push, libmiradac), Cint, (Handle, Handle), l, v)

for f in (:mdac_sdalist_compose_c, :mdac_csdalist_compose, :mdac_csdalist_compose_s)
    @eval $f(m, v, out) = ccall(($(QuoteNode(f)), libmiradac), Cint,
                                (Handle, Handle, Ptr{Handle}), m, v, out)
end

# ---- Multi-env (T7.1) -----------------------------------------------------------------------

for f in (:mdac_nda_import, :mdac_cnda_import, :mdac_sda_import, :mdac_csda_import,
          :mdac_nda_promote_to)
    @eval $f(e, v, out) = ccall(($(QuoteNode(f)), libmiradac), Cint, (Handle, Handle, Ptr{Handle}),
                                e, v, out)
end

# ---- Deferred frees (A.5 item 2) ------------------------------------------------------------
#
# A finalizer only queues the pointer. A pool has no lock and one env is used by one thread at a
# time, so a drain frees only the vectors of the env the calling operation uses (and those of
# retired envs, whose free touches no pool); the others are parked by env until a later drain
# for their env (the next one after a GC, or a retry). `drain!()` frees everything and is for a
# point where no other thread uses an env.
#
# A drain holds the finalizers' spin lock only to swap `ptrs` for an empty vector allocated
# beforehand: no allocation, so no GC, happens under it. Otherwise finalizers run by a GC on
# another thread would find it held and all be put off to a later GC. Neither lock yields: a
# task must not move to another thread in the middle of an operation, since the current env is
# thread-local.

mutable struct FreeQueue{F,E}
    const lock::Threads.SpinLock
    ptrs::Vector{Handle}                            # filled by finalizers
    const parked_lock::Threads.SpinLock             # drains only
    const parked::Dict{Handle,Vector{Handle}}       # env => pointers
    const free::F                                   # the mdac_*_free function
    const envof::E                                  # the mdac_*_env function
end

FreeQueue(free, envof) = FreeQueue(Threads.SpinLock(), Handle[], Threads.SpinLock(),
                                   Dict{Handle,Vector{Handle}}(), free, envof)

const NDA_QUEUE = FreeQueue(mdac_nda_free, mdac_nda_env)
const LIST_QUEUE = FreeQueue(mdac_ndalist_free, mdac_ndalist_env)   # an empty list has no env
const CNDA_QUEUE = FreeQueue(mdac_cnda_free, mdac_cnda_env)
const CLIST_QUEUE = FreeQueue(mdac_cndalist_free, mdac_cndalist_env)
const EXPR_QUEUE = FreeQueue(mdac_expr_free, Returns(C_NULL))       # no env: always freed
const SDA_QUEUE = FreeQueue(mdac_sda_free, mdac_sda_env)
const SLIST_QUEUE = FreeQueue(mdac_sdalist_free, mdac_sdalist_env)
const CSDA_QUEUE = FreeQueue(mdac_csda_free, mdac_csda_env)
const CSLIST_QUEUE = FreeQueue(mdac_csdalist_free, mdac_csdalist_env)
const FREE_QUEUES = (NDA_QUEUE, LIST_QUEUE, CNDA_QUEUE, CLIST_QUEUE, EXPR_QUEUE, SDA_QUEUE, SLIST_QUEUE,
                     CSDA_QUEUE, CSLIST_QUEUE)

"""The finalizer of every handle type: queue the pointer, or retry at the next GC."""
function defer_free(q::FreeQueue, obj)
    if trylock(q.lock)
        push!(q.ptrs, obj.ptr)
        unlock(q.lock)
    else
        finalizer(o -> defer_free(q, o), obj)
    end
    return nothing
end

# `x` gives the env of the operation: a DA object or an env handle; `nothing` is the current env.
op_env(x::Handle) = x
op_env(::Nothing) = (r = Ref{Handle}(C_NULL); mdac_env_current(r); r[])

@inline function drain!(q::FreeQueue, x)
    isempty(q.ptrs) && return nothing
    drain_env!(q, op_env(x))
end

function take_queued!(q::FreeQueue)
    fresh = Handle[]
    lock(q.lock)
    ps = q.ptrs
    q.ptrs = fresh
    unlock(q.lock)
    return ps
end

@noinline function drain_env!(q::FreeQueue, e::Handle)
    ps = take_queued!(q)
    mine = Handle[]
    @lock q.parked_lock begin
        for p in ps
            push!(get!(Vector{Handle}, q.parked, q.envof(p)), p)
        end
        for pe in collect(keys(q.parked))
            (pe == e || pe == C_NULL || mdac_env_retired(pe) != 0) && append!(mine, pop!(q.parked, pe))
        end
    end
    foreach(q.free, mine)
    return nothing
end

"""
    MiraDAC.drain!()

Free every queued vector of every env now. Call it only where no other thread uses an env.
"""
function drain!()
    for q in FREE_QUEUES
        foreach(q.free, take_queued!(q))
        @lock q.parked_lock begin
            foreach(ps -> foreach(q.free, ps), values(q.parked))
            empty!(q.parked)
        end
    end
end

function drain_all!(x)
    e = op_env(x)
    foreach(q -> drain_env!(q, e), FREE_QUEUES)
end

# ---- Allocation with retry on pool exhaustion (A.5 item 3) ----------------------------------

"""Number of times an allocating call retried after `MDAC_ERR_POOL`."""
const POOL_RETRIES = Threads.Atomic{Int}(0)

"""
    alloc_call(f, x=nothing, d=nothing)

Run the allocating C call `f() -> status` for an operation in the env of `x` (see `op_env`):
drain the free queues first; on `MDAC_ERR_POOL` collect garbage and retry, first with an
incremental then with full collections (`retry_rounds()` in all), before throwing
`PoolExhaustedError`. Other failures
throw as `check(st, d)` does.
"""
@inline function alloc_call(f, x=nothing, d=nothing)
    drain!(NDA_QUEUE, x)
    st = f()
    st == MDAC_OK && return nothing
    st == MDAC_ERR_POOL || throw_status(st, d)
    return retry_alloc(f, x, d)
end

# With several threads, `GC.gc` returns without collecting when another thread is already
# collecting (maybe only the young generation), and that thread runs the finalizers: so the
# retry repeats full collections and gives those finalizers a moment. `systemsleep` does not
# yield, so the task stays on this thread (and its current env).
retry_rounds() = Threads.nthreads() == 1 ? 2 : 12

@noinline function retry_alloc(f, x, d)
    for k in 1:retry_rounds()
        Threads.atomic_add!(POOL_RETRIES, 1)
        drain_all!(x)
        GC.gc(k > 1)
        k > 2 && Libc.systemsleep(0.001)
        drain_all!(x)
        st = f()
        st == MDAC_OK && return nothing
        st == MDAC_ERR_POOL || throw_status(st, d)
    end
    throw(PoolExhaustedError(mdac_last_error()))
end
