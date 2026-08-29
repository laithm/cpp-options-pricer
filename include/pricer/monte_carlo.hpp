#pragma once
#include "pricer/payoff.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
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
    detail::validate_option_spec(o);
    if (n_pairs <= 1)
        throw std::invalid_argument("Monte Carlo antithetic pair count must exceed one");
    if (o.T == 0.0)
        return McResult{payoff(o.type, o.S, o.K), 0.0};

    std::mt19937_64 gen(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    double drift = (o.r - o.q - 0.5 * o.sigma * o.sigma) * o.T;
    double diffusion = o.sigma * std::sqrt(o.T);
    double disc = std::exp(-o.r * o.T);
    if (!std::isfinite(drift) || !std::isfinite(diffusion) || !std::isfinite(disc))
        throw std::invalid_argument("Monte Carlo inputs produce invalid GBM factors");

    if (o.sigma == 0.0) {
        double terminal_spot = o.S * std::exp(drift);
        if (!std::isfinite(terminal_spot))
            throw std::invalid_argument(
                "Monte Carlo inputs produce a non-finite deterministic value");
        double price = disc * payoff(o.type, terminal_spot, o.K);
        if (!std::isfinite(price))
            throw std::invalid_argument(
                "Monte Carlo inputs produce a non-finite deterministic value");
        return McResult{price, 0.0};
    }

    double sum = 0.0, sum_sq = 0.0;
    for (long i = 0; i < n_pairs; ++i) {
        double z = Z(gen);
        double st_up = o.S * std::exp(drift + diffusion * z);
        double st_dn = o.S * std::exp(drift - diffusion * z);
        if (!std::isfinite(st_up) || !std::isfinite(st_dn))
            throw std::invalid_argument(
                "Monte Carlo inputs produce non-finite terminal values");
        double pay = 0.5 * (payoff(o.type, st_up, o.K) + payoff(o.type, st_dn, o.K));
        double sample = disc * pay;
        if (!std::isfinite(sample))
            throw std::invalid_argument("Monte Carlo inputs produce non-finite payoffs");
        sum += sample;
        sum_sq += sample * sample;
    }
    double n = static_cast<double>(n_pairs);
    double mean = sum / n;
    double second_moment = sum_sq / n;
    double variance_numerator = second_moment - mean * mean;
    double roundoff_tolerance = 16.0 * std::numeric_limits<double>::epsilon()
                              * std::max(second_moment, mean * mean);
    if (variance_numerator < 0.0
        && -variance_numerator <= roundoff_tolerance)
        variance_numerator = 0.0;
    double var = variance_numerator * (n / (n - 1.0)); // sample variance
    if (!std::isfinite(mean) || !std::isfinite(var) || var < 0.0)
        throw std::invalid_argument("Monte Carlo inputs produce an invalid variance estimate");
    return McResult{mean, std::sqrt(var / n)};
}

} // namespace pricer
