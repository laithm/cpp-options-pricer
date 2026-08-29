#include "pricer/black_scholes.hpp"
#include "pricer/binomial.hpp"
#include "pricer/monte_carlo.hpp"
#include <cstdio>
#include <fstream>
#include <utility>
#include <vector>

using namespace pricer;

int main() {
    OptionSpec o{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    double bs = bs_price(o);

    std::vector<int> crr_steps{10, 50, 100, 500, 1000};
    std::vector<long> mc_pairs{10000, 100000, 1000000};

    std::printf("Reference: S=%.0f K=%.0f r=%.2f q=%.2f sigma=%.2f T=%.1f (call)\n",
                o.S, o.K, o.r, o.q, o.sigma, o.T);
    std::printf("Black-Scholes: %.6f\n\n", bs);

    std::printf("%-8s %-10s %-12s %-12s\n", "method", "param", "price", "abs_err");
    for (int n : crr_steps)
        std::printf("%-8s %-10d %-12.6f %-12.2e\n", "CRR", n,
                    crr_price(o, n, false), std::fabs(crr_price(o, n, false) - bs));
    std::vector<McResult> mc_results;
    for (long n : mc_pairs) {
        McResult m = mc_price(o, n, 42);
        mc_results.push_back(m);
        std::printf("%-8s %-10ld %-12.6f se=%.4f\n", "MC", n, m.price, m.std_error);
    }

    std::printf("\nGreek    analytic     finite_diff\n");
    auto fd = [&](char w, double h) {
        OptionSpec a = o, b = o;
        switch (w) { case 'S': a.S+=h; b.S-=h; break; case 'v': a.sigma+=h; b.sigma-=h; break;
                     case 'r': a.r+=h; b.r-=h; break; case 't': a.T+=h; b.T-=h; break; }
        return std::make_pair(bs_price(a), bs_price(b));
    };
    double h = 1e-4;
    struct G { const char* name; double an; double fdv; };
    auto p_s = fd('S', h); auto p_v = fd('v', h); auto p_r = fd('r', h); auto p_t = fd('t', h);
    std::vector<G> greeks{
        {"delta", bs_delta(o), (p_s.first - p_s.second) / (2*h)},
        {"gamma", bs_gamma(o), (p_s.first - 2*bs + p_s.second) / (h*h)},
        {"vega",  bs_vega(o),  (p_v.first - p_v.second) / (2*h)},
        {"rho",   bs_rho(o),   (p_r.first - p_r.second) / (2*h)},
        {"theta", bs_theta(o), -(p_t.first - p_t.second) / (2*h)},
    };
    for (auto& g : greeks)
        std::printf("%-8s %-12.6f %-12.6f\n", g.name, g.an, g.fdv);

    // Write results.json (hand-transcribed into the site later).
    std::ofstream f("output/results.json");
    f << "{\n";
    f << "  \"contract\": {\"S\":100,\"K\":100,\"r\":0.05,\"q\":0.0,\"sigma\":0.20,\"T\":1.0,\"type\":\"call\"},\n";
    f << "  \"bs\": " << bs << ",\n";
    f << "  \"convergence\": [\n";
    bool first = true;
    for (int n : crr_steps) {
        f << (first ? "    " : ",\n    ") << "{\"method\":\"CRR\",\"param\":" << n
          << ",\"price\":" << crr_price(o, n, false) << "}";
        first = false;
    }
    for (size_t i = 0; i < mc_pairs.size(); ++i) {
        f << ",\n    {\"method\":\"MC\",\"param\":" << mc_pairs[i]
          << ",\"price\":" << mc_results[i].price
          << ",\"std_error\":" << mc_results[i].std_error << "}";
    }
    f << "\n  ],\n  \"greeks\": [\n";
    for (size_t i = 0; i < greeks.size(); ++i)
        f << (i ? ",\n    " : "    ") << "{\"name\":\"" << greeks[i].name
          << "\",\"analytic\":" << greeks[i].an
          << ",\"finite_diff\":" << greeks[i].fdv << "}";
    f << "\n  ]\n}\n";
    f.close();
    std::printf("\nwrote output/results.json\n");
    return 0;
}
