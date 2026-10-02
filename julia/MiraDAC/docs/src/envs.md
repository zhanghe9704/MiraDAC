# Multiple envs

`init!(order, nvars, poolsize)` makes the *default* env and makes it current; `clear!()`
retires it. `DAEnv(order, nvars, poolsize; table=false)` makes another env without changing the
current one.

```julia
init!(4, 3, 10_000)
e = DAEnv(10, 2, 500)
q = with_env(e) do                  # e is current inside, restored afterwards
    davar(1) * davar(2)
end
p = davar(1; env=e)                 # or name the env with the `env` keyword
current_env() == default_env()      # true
e.order, e.nvars, e.count           # properties
```

- Every operation runs in the env of its first DA operand, whichever env is current. Vectors of
  different envs in one operation throw `EnvError`.
- `import_vec(env, v)` copies `v` into `env` (`ArgumentError` if the layouts differ);
  `promote_sda(env, v)` promotes an `NDA` or `CNDA` into `env` as an `SDA` or `CSDA`.
- `close(env)` releases an env. Using a vector of a closed (or cleared) env throws `EnvError`;
  dropping it is safe, its memory is freed when it is collected.
- `with_order(f, env, n)` truncates one env only.
