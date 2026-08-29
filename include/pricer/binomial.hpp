#pragma once
#include "pricer/payoff.hpp"
#include <cmath>
#include <vector>

namespace pricer {

// Cox-Ross-Rubinstein binomial price. american=true allows early exercise.
inline double crr_price(const OptionSpec& o, int n, bool american) {
    double dt = o.T / n;
    double u = std::exp(o.sigma * std::sqrt(dt));
    double d = 1.0 / u;
    double disc = std::exp(-o.r * dt);
    double p = (std::exp((o.r - o.q) * dt) - d) / (u - d);

    // Terminal payoffs: node j has j up-moves, (n-j) down-moves.
    std::vector<double> v(n + 1);
    for (int j = 0; j <= n; ++j) {
        double S_T = o.S * std::pow(u, j) * std::pow(d, n - j);
        v[j] = payoff(o.type, S_T, o.K);
    }

    // Backward induction.
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            double cont = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (american) {
                double S_ij = o.S * std::pow(u, j) * std::pow(d, i - j);
                cont = std::max(cont, payoff(o.type, S_ij, o.K));
            }
            v[j] = cont;
        }
    }
    return v[0];
}

} // namespace pricer
