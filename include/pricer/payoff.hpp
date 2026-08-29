#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pricer {

enum class OptionType { Call, Put };

struct OptionSpec {
    double S;      // spot
    double K;      // strike
    double r;      // risk-free rate (cont. comp.)
    double q;      // continuous dividend yield
    double sigma;  // volatility
    double T;      // time to maturity (years)
    OptionType type;
};

namespace detail {

inline void validate_option_spec(const OptionSpec& o) {
    if (!std::isfinite(o.S) || o.S <= 0.0)
        throw std::invalid_argument("spot must be finite and positive");
    if (!std::isfinite(o.K) || o.K <= 0.0)
        throw std::invalid_argument("strike must be finite and positive");
    if (!std::isfinite(o.r))
        throw std::invalid_argument("risk-free rate must be finite");
    if (!std::isfinite(o.q))
        throw std::invalid_argument("dividend yield must be finite");
    if (!std::isfinite(o.sigma) || o.sigma < 0.0)
        throw std::invalid_argument("volatility must be finite and non-negative");
    if (!std::isfinite(o.T) || o.T < 0.0)
        throw std::invalid_argument("maturity must be finite and non-negative");
    if (o.type != OptionType::Call && o.type != OptionType::Put)
        throw std::invalid_argument("option type is invalid");
}

} // namespace detail

// Intrinsic payoff at a terminal spot S_T.
inline double payoff(OptionType type, double S_T, double K) {
    return type == OptionType::Call ? std::max(S_T - K, 0.0)
                                    : std::max(K - S_T, 0.0);
}

} // namespace pricer
