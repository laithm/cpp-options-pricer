#pragma once
#include "pricer/payoff.hpp"
#include <cstddef>
#include <cmath>
#include <vector>

namespace pricer {

// Cox-Ross-Rubinstein binomial price. american=true allows early exercise.
inline double crr_price(const OptionSpec& o, int n, bool american) {
    detail::validate_option_spec(o);
    if (n <= 0)
        throw std::invalid_argument("CRR steps must be positive");
    if (o.T <= 0.0)
        throw std::invalid_argument("CRR requires positive maturity");
    if (o.sigma <= 0.0)
        throw std::invalid_argument("CRR requires positive volatility");

    double dt = o.T / static_cast<double>(n);
    double u = std::exp(o.sigma * std::sqrt(dt));
    double d = 1.0 / u;
    double disc = std::exp(-o.r * dt);
    double growth = std::exp((o.r - o.q) * dt);
    double denominator = u - d;
    if (!std::isfinite(dt) || !std::isfinite(u) || !std::isfinite(d)
        || !std::isfinite(disc) || !std::isfinite(growth)
        || !std::isfinite(denominator) || denominator <= 0.0)
        throw std::invalid_argument("CRR inputs produce invalid tree factors");

    double p = (growth - d) / denominator;
    if (!std::isfinite(p) || p < 0.0 || p > 1.0)
        throw std::invalid_argument(
            "CRR risk-neutral probability must be finite and in [0, 1]");

    // Terminal payoffs: node j has j up-moves, (n-j) down-moves.
    std::vector<double> v(static_cast<std::size_t>(n) + 1U);
    for (int j = 0; j <= n; ++j) {
        double S_T = o.S * std::pow(u, j) * std::pow(d, n - j);
        if (!std::isfinite(S_T))
            throw std::invalid_argument("CRR inputs produce non-finite node values");
        v[j] = payoff(o.type, S_T, o.K);
    }

    // Backward induction.
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            double cont = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (american) {
                double S_ij = o.S * std::pow(u, j) * std::pow(d, i - j);
                if (!std::isfinite(S_ij))
                    throw std::invalid_argument(
                        "CRR inputs produce non-finite node values");
                cont = std::max(cont, payoff(o.type, S_ij, o.K));
            }
            v[j] = cont;
        }
    }
    return v[0];
}

} // namespace pricer
