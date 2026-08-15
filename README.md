# Numerical Ritz Method (1D)

A small Julia package implementing the Ritz method for 1D calculus-of-variations problems. It was built as a personal learning project to accompany Chapter 8 of Gelfand & Fomin's *Calculus of Variations*, mainly to numerically check coefficients I'd worked out by hand.

## Background

Given a functional

```math
J[y] = \int_{a}^{b} F(x, y(x), y'(x))\, dx
```

the Ritz method approximates the minimizing function $y(x)$ by restricting the search to a trial solution built from a finite set of basis functions,

```math
y(x) = \phi_0(x) + \sum_{i=1}^{n} a_i \phi_i(x)
```

where $\phi_0$ satisfies any inhomogeneous boundary conditions and each $\phi_i$ vanishes at the boundary. This turns the variational problem into a finite-dimensional optimization over the coefficients $a_i$, which this package solves by evaluating $J[y]$ via Gauss-Legendre quadrature and minimizing it with [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/).

## Installation

The package isn't registered, so install it directly from this repo:

```julia
using Pkg
Pkg.develop(url="https://github.com/shankal17/Numerical-Ritz-Method-1D.git")
```

Or, working from a local clone:

```julia
using Pkg
Pkg.develop(path="path/to/Numerical-Ritz-Method-1D")
```

## API

- **`Interval(low, high)`** — a closed domain `[low, high]`.

- **`minimizefunctional(F, domain, phi_s; phi_0, num_quad_nodes, initial_coeffs, method, options)`** — approximates the minimizer of `J[y]` over `domain`, given:
  - `F` — the integrand, called as `F(x, y, y')`
  - `domain::Interval` — `[a, b]`
  - `phi_s::Vector` — basis functions `φᵢ` (callable)
  - `phi_0` (kwarg, default `x -> 0.0`) — fixed term for boundary conditions
  - `num_quad_nodes` (kwarg, default `15`) — Number of Gauss-Legendre quadrature nodes
  - `initial_coeffs` (kwarg, default `zeros(length(phi_s))`)
  - `method` (kwarg, default `Optim.BFGS()`)
  - `options` (kwarg, default `Optim.Options()`)

  Returns a `RitzResult`.

- **`RitzResult`** — fields `coeffs`, `minimized_value`, `domain`, `converged`, `y` (the trial solution as a callable, e.g. `result.y(0.5)`).

- **`summarize(result)`** — human-readable summary string (also used by `Base.show`, so `println(result)` works directly).

## Usage

```julia
using RitzMethod
import Optim

# J[y] = ∫ (y'^2 - y^2 - 2xy) dx over [0, 1]
integrand(x, y, y_prime) = y_prime^2 - y^2 - 2 * x * y
domain = Interval(0.0, 1.0)

# Basis functions vanishing at both endpoints
basis = [x -> x^k * (1 - x) for k in 1:3]

result = minimizefunctional(integrand, domain, basis; num_quad_nodes=5)

println(summarize(result))
result.coeffs      # optimized coefficients c_i
result.y(0.5)       # approximate solution at x = 0.5
```

## Worked examples

`scripts/` contains full worked examples for two Gelfand & Fomin problems, each comparing 1-, 2-, and 3-term Ritz approximations against the closed-form exact solution:

- [`scripts/8.2/report.jmd`](scripts/8.2/report.jmd) — Problem 8.2, domain `[0, 1]`
- [`scripts/8.4/report.jmd`](scripts/8.4/report.jmd) — Problem 8.4, domain `[0, 2]`

These are [Weave.jl](https://weavejl.mpastell.com/stable/) documents. Run `scripts/build_reports.jl` to render them to HTML (output goes to `scripts/*/report/report.html`).

## Project structure

```
src/
  RitzMethod.jl     # the package: Interval, minimizefunctional, RitzResult, summarize
scripts/             # separate Julia environment with worked examples
  8.2/report.jmd
  8.4/report.jmd
  build_reports.jl  # renders all scripts/*/report.jmd to HTML
```

## Failure Modes/Gottchas
There are a few points that the user should be aware of

- **Nonconvex Objective Functions** — When using the `BFGS` optimizer (the default), nonconvex objective functions can either converge on a local minimum (will still report `converged = true`), or never converge.
- **No boundary-condition enforcement** — The code does not check if the basis functions actually satisfy the boundary conditions. The user must make sure that their chosen basis functions make sense.
- **Ill-conditioning from basis choice** — As more terms are added, the optimization landscape might become ill-conditioned, degrading BFGS's Hessian approximation.
- **Fixed-order quadrature** — The user should make sure that their choice of `num_quad_nodes` is appropriate.
- **First-derivative ceiling** — only integrands of form `F(x, y, y')` are supported.
- **Automatic-differentiation fragility** — both the basis-function derivatives and the objective gradient rely on `ForwardDiff`; any `F` or `phi_i` using non-AD-friendly operations/logic will error or produce incorrect gradients.
