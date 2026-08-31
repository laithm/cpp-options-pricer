# C++ Options Pricer

[![CI](https://github.com/laithm/cpp-options-pricer/actions/workflows/ci.yml/badge.svg)](https://github.com/laithm/cpp-options-pricer/actions/workflows/ci.yml)

I built this small C++20 pricer to make several pricing methods check each
other on the same contract. It supports European calls and puts,
early-exercise handling in the binomial tree, and Monte Carlo estimates with a
reported standard error. The point is to understand the numerical methods and
their failure checks, not to pretend this is a production pricing system.

## Methods

- **Black-Scholes** — closed-form European call and put prices with continuous dividend yield.
- **Cox-Ross-Rubinstein (CRR)** — recombining binomial tree for European and American exercise.
- **Monte Carlo** — terminal-value simulation under risk-neutral geometric Brownian motion using antithetic variates. The estimator returns the discounted price and its standard error.

## Greeks

Black-Scholes Delta, Gamma, Vega, Rho and Theta are implemented analytically. Vega and Rho are reported per unit change in volatility and interest rate; Theta is annual calendar decay.

The validation harness and tests compare the analytic values with centered finite differences. Gamma uses a centered second derivative, while Theta is checked against the negative derivative with respect to time to maturity.

## Validation

The automated tests use the reference contract \(S=K=100\), \(r=5\%\), \(q=0\), \(\sigma=20\%\) and \(T=1\).

| Check | Criterion |
| --- | --- |
| Textbook Black-Scholes values | Call `10.450583` and put `5.573526`, each within `1e-4` |
| Put-call parity | Residual within `1e-8` |
| CRR convergence | 1,000-step European call and put prices within `5e-2` of Black-Scholes |
| Exercise invariant | American put value is not below the corresponding European put value |
| Monte Carlo error | 200,000 antithetic pairs; estimate within three reported standard errors of Black-Scholes |
| Greeks | Analytic values agree with centered finite differences under explicit per-Greek tolerances |

The executable also reports a CRR step sequence, seeded Monte Carlo estimates and analytic-versus-finite-difference Greeks, then exports the results as JSON.

## Build

The commands below require CMake 3.20 or later and a C++20 compiler.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Run the validation harness

From the repository root:

```bash
mkdir -p output
./build/pricer_run
```

This writes `output/results.json`. Generated build and output files are not committed.

## Layout

- `include/pricer/payoff.hpp` — option specification and terminal payoff
- `include/pricer/black_scholes.hpp` — Black-Scholes prices and analytic Greeks
- `include/pricer/binomial.hpp` — European and American CRR valuation
- `include/pricer/monte_carlo.hpp` — antithetic GBM Monte Carlo and standard error
- `src/main.cpp` — validation and JSON export harness
- `tests/test_pricer.cpp` — numerical consistency and invariant checks

## Design notes

- Pricing routines are implemented in headers; the executable is a validation and JSON-export harness.
- Monte Carlo accepts an explicit seed. Tests and the validation harness use fixed seeds for repeatable runs within a given standard-library implementation.
- Closed-form, tree, simulation and finite-difference calculations are cross-validated against one another.
- Monte Carlo estimates include an explicit standard error.

## Limitations

- Constant volatility
- GBM/lognormal dynamics
- Deterministic interest rates and dividend yield
- No volatility-surface calibration
- No market-data ingestion
- No execution functionality
- Educational and research implementation, not production pricing infrastructure

## What I would add next

- Implied-volatility solver
- A convergence and runtime study across tree depth and Monte Carlo sample size
