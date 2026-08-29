#pragma once
#include <algorithm>

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

// Intrinsic payoff at a terminal spot S_T.
inline double payoff(OptionType type, double S_T, double K) {
    return type == OptionType::Call ? std::max(S_T - K, 0.0)
                                    : std::max(K - S_T, 0.0);
}

} // namespace pricer
