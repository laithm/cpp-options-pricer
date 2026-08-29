#include "pricer/payoff.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/binomial.hpp"
#include "pricer/monte_carlo.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace pricer;

static int failures = 0;
#define CHECK_NEAR(a, b, tol, name)                                      \
    do {                                                                 \
        double _a = (a), _b = (b);                                       \
        if (std::fabs(_a - _b) > (tol)) {                                \
            std::printf("FAIL %s: %.6f vs %.6f (tol %.2g)\n",            \
                        name, _a, _b, (double)(tol));                    \
            ++failures;                                                  \
        } else {                                                         \
            std::printf("ok   %s\n", name);                              \
        }                                                                \
    } while (0)

static void test_payoff() {
    CHECK_NEAR(payoff(OptionType::Call, 110, 100), 10.0, 1e-12, "call payoff");
    CHECK_NEAR(payoff(OptionType::Put, 90, 100), 10.0, 1e-12, "put payoff");
    CHECK_NEAR(payoff(OptionType::Call, 90, 100), 0.0, 1e-12, "call OTM payoff");
}

static void test_black_scholes() {
    // S=100,K=100,r=0.05,q=0,sigma=0.20,T=1 -> textbook call 10.4506, put 5.5735
    OptionSpec call{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    OptionSpec put = call; put.type = OptionType::Put;
    CHECK_NEAR(bs_price(call), 10.450583, 1e-4, "bs call price");
    CHECK_NEAR(bs_price(put), 5.573526, 1e-4, "bs put price");

    // Put-call parity: C - P = S e^{-qT} - K e^{-rT}
    double parity = call.S * std::exp(-call.q * call.T)
                  - call.K * std::exp(-call.r * call.T);
    CHECK_NEAR(bs_price(call) - bs_price(put), parity, 1e-8, "put-call parity");
}

static double bump_price(OptionSpec o, char which, double h) {
    switch (which) {
        case 'S': o.S += h; break;
        case 'v': o.sigma += h; break;
        case 't': o.T += h; break;   // forward bump in time
        case 'r': o.r += h; break;
    }
    return bs_price(o);
}

static void test_greeks() {
    OptionSpec o{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    double h = 1e-4;
    double fd_delta = (bump_price(o, 'S', h) - bump_price(o, 'S', -h)) / (2 * h);
    double fd_gamma = (bump_price(o, 'S', h) - 2 * bs_price(o) + bump_price(o, 'S', -h)) / (h * h);
    double fd_vega  = (bump_price(o, 'v', h) - bump_price(o, 'v', -h)) / (2 * h);
    double fd_rho   = (bump_price(o, 'r', h) - bump_price(o, 'r', -h)) / (2 * h);
    // theta = -d(price)/dt ; finite-diff via time bump
    double fd_theta = -(bump_price(o, 't', h) - bump_price(o, 't', -h)) / (2 * h);

    CHECK_NEAR(bs_delta(o), fd_delta, 1e-4, "delta vs FD");
    CHECK_NEAR(bs_gamma(o), fd_gamma, 1e-2, "gamma vs FD");
    CHECK_NEAR(bs_vega(o),  fd_vega,  1e-3, "vega vs FD");
    CHECK_NEAR(bs_rho(o),   fd_rho,   1e-3, "rho vs FD");
    CHECK_NEAR(bs_theta(o), fd_theta, 1e-3, "theta vs FD");
}

static void test_binomial() {
    OptionSpec call{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    // European CRR converges to BS as steps grow
    CHECK_NEAR(crr_price(call, 1000, false), bs_price(call), 5e-2, "CRR->BS call");

    OptionSpec put = call; put.type = OptionType::Put;
    CHECK_NEAR(crr_price(put, 1000, false), bs_price(put), 5e-2, "CRR->BS put");

    // American put >= European put (early exercise premium, never negative)
    double amer = crr_price(put, 1000, true);
    double euro = crr_price(put, 1000, false);
    if (amer < euro - 1e-6) { std::printf("FAIL amer put < euro put\n"); ++failures; }
    else std::printf("ok   american put >= european put\n");
}

static void test_monte_carlo() {
    OptionSpec call{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    McResult mc = mc_price(call, 200000, 12345);
    double bs = bs_price(call);
    // Price within 3 standard errors of BS (honest error bar).
    if (std::fabs(mc.price - bs) > 3 * mc.std_error) {
        std::printf("FAIL MC within 3 SE: |%.4f-%.4f|=%.4f > 3*%.4f\n",
                    mc.price, bs, std::fabs(mc.price - bs), mc.std_error);
        ++failures;
    } else {
        std::printf("ok   MC within 3 SE of BS (price %.4f, se %.4f)\n",
                    mc.price, mc.std_error);
    }
    if (mc.std_error <= 0) { std::printf("FAIL MC std_error <= 0\n"); ++failures; }
}

int main() {
    test_payoff();
    test_black_scholes();
    test_greeks();
    test_binomial();
    test_monte_carlo();
    if (failures) { std::printf("\n%d FAILURES\n", failures); return 1; }
    std::printf("\nALL PASS\n");
    return 0;
}
