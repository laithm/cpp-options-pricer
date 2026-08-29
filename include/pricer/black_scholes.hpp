#pragma once
#include "pricer/payoff.hpp"
#include <cmath>

namespace pricer {

inline double norm_pdf(double x) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

namespace detail {

inline void validate_black_scholes_inputs(const OptionSpec& o) {
    validate_option_spec(o);
    if (o.T <= 0.0)
        throw std::invalid_argument("Black-Scholes requires positive maturity");
    if (o.sigma <= 0.0)
        throw std::invalid_argument("Black-Scholes requires positive volatility");
}

inline double bs_d1_unchecked(const OptionSpec& o) {
    return (std::log(o.S / o.K) + (o.r - o.q + 0.5 * o.sigma * o.sigma) * o.T)
           / (o.sigma * std::sqrt(o.T));
}

inline double bs_d2_unchecked(const OptionSpec& o) {
    return bs_d1_unchecked(o) - o.sigma * std::sqrt(o.T);
}

} // namespace detail

inline double d1(const OptionSpec& o) {
    detail::validate_black_scholes_inputs(o);
    return detail::bs_d1_unchecked(o);
}

inline double d2(const OptionSpec& o) {
    detail::validate_black_scholes_inputs(o);
    return detail::bs_d2_unchecked(o);
}

inline double bs_price(const OptionSpec& o) {
    detail::validate_black_scholes_inputs(o);
    double D1 = detail::bs_d1_unchecked(o), D2 = detail::bs_d2_unchecked(o);
    double disc_q = std::exp(-o.q * o.T);
    double disc_r = std::exp(-o.r * o.T);
    if (o.type == OptionType::Call)
        return o.S * disc_q * norm_cdf(D1) - o.K * disc_r * norm_cdf(D2);
    return o.K * disc_r * norm_cdf(-D2) - o.S * disc_q * norm_cdf(-D1);
}

inline double bs_delta(const OptionSpec& o) {
    detail::validate_black_scholes_inputs(o);
    double disc_q = std::exp(-o.q * o.T);
    double D1 = detail::bs_d1_unchecked(o);
    return o.type == OptionType::Call ? disc_q * norm_cdf(D1)
                                      : disc_q * (norm_cdf(D1) - 1.0);
}
inline double bs_gamma(const OptionSpec& o) {
    detail::validate_black_scholes_inputs(o);
    return std::exp(-o.q * o.T) * norm_pdf(detail::bs_d1_unchecked(o))
           / (o.S * o.sigma * std::sqrt(o.T));
}
inline double bs_vega(const OptionSpec& o) {  // per 1.0 change in sigma
    detail::validate_black_scholes_inputs(o);
    return o.S * std::exp(-o.q * o.T) * norm_pdf(detail::bs_d1_unchecked(o))
           * std::sqrt(o.T);
}
inline double bs_rho(const OptionSpec& o) {   // per 1.0 change in r
    detail::validate_black_scholes_inputs(o);
    double disc_r = std::exp(-o.r * o.T);
    double D2 = detail::bs_d2_unchecked(o);
    return o.type == OptionType::Call
        ? o.K * o.T * disc_r * norm_cdf(D2)
        : -o.K * o.T * disc_r * norm_cdf(-D2);
}
inline double bs_theta(const OptionSpec& o) { // per 1.0 year (calendar decay)
    detail::validate_black_scholes_inputs(o);
    double D1 = detail::bs_d1_unchecked(o), D2 = detail::bs_d2_unchecked(o);
    double disc_q = std::exp(-o.q * o.T), disc_r = std::exp(-o.r * o.T);
    double term1 = -o.S * disc_q * norm_pdf(D1) * o.sigma / (2 * std::sqrt(o.T));
    if (o.type == OptionType::Call)
        return term1 - o.r * o.K * disc_r * norm_cdf(D2)
               + o.q * o.S * disc_q * norm_cdf(D1);
    return term1 + o.r * o.K * disc_r * norm_cdf(-D2)
           - o.q * o.S * disc_q * norm_cdf(-D1);
}

} // namespace pricer
