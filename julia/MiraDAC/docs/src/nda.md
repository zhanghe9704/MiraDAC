# NDA

An `NDA` is a numeric DA vector of the current env (or of the `env` keyword's).

```julia
init!(4, 3, 10_000)
x = NDA(2.0)                        # a constant
v = NDA([1.0, 2.0, 0.5])            # coefficients in monomial order (see exponents)
x1, x2, x3 = davar(1), davar(2), davar(3)
f = sin(x1) * exp(x2 + x3) / (1 + x1^2)
```

## Arithmetic and functions

`+ - * /` combine `NDA`s and `Real`s on either side; `^` takes an `Integer` or a `Real`
exponent. The math functions are `Base`'s methods: `exp log sqrt sin cos tan asin acos atan sinh
cosh tanh asinh acosh atanh`, plus `MiraDAC.erf` (exported; no SpecialFunctions dependency).
`abs(v)` and `norm(v)` are the largest coefficient magnitude, as in C++. A function outside its
domain at the constant part throws `DomainError`, for example `asin(NDA(2.0))`.

Each operation above returns a new vector; [Performance](@ref) describes the in-place forms
`add!`, `sub!`, `mul!`, `div!` and `exp!`, `sin!`, … .

## Inspection

| | |
|---|---|
| `con(v)` | the constant term |
| `coeff(v, exps)` | the coefficient of the monomial with exponents `exps` |
| `coeffs(v)` | the coefficients in monomial order, up to the last non-zero one |
| `nterms(v)` | the number of non-zero coefficients |
| `exponents(env)` | an `nvars × full_length` matrix: column `k` is the monomial of coefficient `k` |
| `iszero(v; eps)` | all coefficients below `eps` |
| `show` | compact in a container; the full C++ table in the REPL (`display`) |

## Truncation order

`with_order(f, n)` runs `f` with the current env truncated at order `n`, then restores the
previous order (also when `f` throws). Calls nest:

```julia
with_order(2) do
    y = exp(x1)                     # terms above order 2 are dropped
end
```

`get_eps()` and `set_eps!(x)` read and set the threshold below which coefficients are dropped.

## Maps and algorithms

An `NDAList` holds a map in C++ (no per-element copy on the C side); every map argument also
takes a `Vector{NDA}`.

| | |
|---|---|
| `der(v, i)`, `integ(v, i)` | derivative and integral with respect to variable `i` |
| `substitute(v, i, x)` | `x` (a number or an `NDA`) for variable `i`; also several at once and for a map |
| `compose(m, v)` | `m` composed with the vectors `v` (an `NDAList`), or evaluated at a real or complex point |
| `inv_map(m)` | the inverse of a map with zero constant parts (`ArgumentError` if singular) |
| `evaluate_map(m, pts)` | `m` at the columns of the `nvars × N` matrix `pts`: a `length(m) × N` matrix |
| `v(pt)` | `v` at the point `pt` |

```julia
m = NDAList([x1 + x2^2, x2 - x1^3, x3])   # a map of all nvars variables
mi = inv_map(m)
compose(m, mi)                      # the identity map, up to the truncation order
evaluate_map(m, rand(3, 1000))      # 3 × 1000
```
