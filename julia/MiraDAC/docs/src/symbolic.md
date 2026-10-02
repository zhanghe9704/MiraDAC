# Symbolic

Symbolic support needs a library built with SymEngine (`MiraDAC.HAS_SYMBOLIC == true`).

## SymExpr

A `SymExpr` is a symbolic scalar (C++ `SymEngine::Expression`); the name avoids Julia's `Expr`.

```julia
a, b = dasymbols("a b")
e = (a + 2b)^2
MiraDAC.expand(e)                   # 4*a*b + a**2 + 4*b**2
diff(e, a)
MiraDAC.subs(e, Dict(a => 1, b => SymExpr("c")))
Float64(MiraDAC.subs(e, Dict(a => 1, b => 2)))   # 25.0
SymExpr("sin(x)/2")                 # parsed
```

An `Integer` stays an exact integer; another `Real` becomes a floating-point number. `subs`,
`expand` and `free_symbols` are not exported (SymEngine.jl exports the same names); use
`MiraDAC.subs` and so on.

## SDA

An `SDA` is a DA vector with `SymExpr` coefficients.

```julia
init!(4, 3, 10_000)
s = 1 + (1 + a) * davar(1) + b * sdavar(2)     # SymExpr * NDA gives an SDA
t = exp(s) * davar(3)                          # SDA ⊕ NDA gives an SDA
con(t), coeffs(t)                              # SymExprs
v = evaluate(t, Dict(a => 0.3, b => 0.2))      # an NDA
promote_sda(davar(1))                          # NDA → SDA
```

Operators combine `SDA` with `SDA`, `NDA`, `SymExpr` and `Real`; the functions are the NDA set
except `asinh acosh atanh`. `simplify`, `MiraDAC.expand` and `MiraDAC.subs` act on every
coefficient. `der`, `integ`, `substitute` and `compose` work on `SDA`s and `SDAList`s.

## CSDA

A `CSDA` is a complex symbolic DA vector `re + im*i` with `SDA` parts; the imaginary unit stays
outside the coefficients, so the parts stay separable.

```julia
z = CSDA(1.5 + a * sdavar(1), 0.5 + b * sdavar(1))
ez = exp(z)
evaluate(ez, Dict(a => 0.3, b => 0.2))         # a CNDA
```

Its functions are those of `CNDA`; `abs(z)` is the modulus as an `SDA`. Keep the constant parts
of function arguments numeric: a symbolic constant part makes the expressions grow very fast.

## SymEngine.jl

With SymEngine.jl loaded, a package extension converts in both directions through strings:
`SymExpr(x::SymEngine.Basic)` and `SymEngine.Basic(x::SymExpr)`; `SDA` constructors and
`evaluate` keys accept `SymEngine.Basic`. The two packages use different SymEngine libraries, so
there is no zero-copy exchange. Use `import SymEngine` (with `using`, the unqualified name
`coeff` is ambiguous).
