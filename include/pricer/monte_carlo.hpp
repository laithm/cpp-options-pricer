#pragma once
#include "pricer/payoff.hpp"
#include <cmath>
#include <random>

namespace pricer {

struct McResult {
    double price;
    double std_error;
};

// Terminal-value GBM Monte Carlo with antithetic variates.
// n_pairs antithetic pairs are drawn; each pair's mean discounted payoff is one
// sample, so the standard error reflects the variance reduction honestly.
inline McResult mc_price(const OptionSpec& o, long n_pairs, unsigned seed) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    double drift = (o.r - o.q - 0.5 * o.sigma * o.sigma) * o.T;
    double diffusion = o.sigma * std::sqrt(o.T);
    double disc = std::exp(-o.r * o.T);

    double sum = 0.0, sum_sq = 0.0;
    for (long i = 0; i < n_pairs; ++i) {
        double z = Z(gen);
        double st_up = o.S * std::exp(drift + diffusion * z);
        double st_dn = o.S * std::exp(drift - diffusion * z);
        double pay = 0.5 * (payoff(o.type, st_up, o.K) + payoff(o.type, st_dn, o.K));
        double sample = disc * pay;
        sum += sample;
        sum_sq += sample * sample;
    }
    double n = static_cast<double>(n_pairs);
    double mean = sum / n;
    double var = (sum_sq / n - mean * mean) * (n / (n - 1.0)); // sample variance
    return McResult{mean, std::sqrt(var / n)};
}

} // namespace pricer
