#include "pricer/payoff.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/binomial.hpp"
#include "pricer/monte_carlo.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

using namespace pricer;

static int failures = 0;
#define CHECK_NEAR(a, b, tol, name)                                      \
    do {                                                                 \
        double _a = (a), _b = (b);                                       \
        if (!std::isfinite(_a) || !std::isfinite(_b)                     \
            || std::fabs(_a - _b) > (tol)) {                             \
            std::printf("FAIL %s: %.6f vs %.6f (tol %.2g)\n",            \
                        name, _a, _b, (double)(tol));                    \
            ++failures;                                                  \
        } else {                                                         \
            std::printf("ok   %s\n", name);                              \
        }                                                                \
    } while (0)

template <typename Fn>
static void check_invalid_argument(Fn&& fn, const char* name) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        std::printf("ok   %s\n", name);
        return;
    } catch (const std::exception& e) {
        std::printf("FAIL %s: wrong exception (%s)\n", name, e.what());
        ++failures;
        return;
    } catch (...) {
        std::printf("FAIL %s: wrong non-standard exception\n", name);
        ++failures;
        return;
    }
    std::printf("FAIL %s: expected std::invalid_argument\n", name);
    ++failures;
}

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
    // Check the estimate against a three-standard-error interval around BS.
    if (!std::isfinite(mc.price) || !std::isfinite(mc.std_error)
        || std::fabs(mc.price - bs) > 3 * mc.std_error) {
        std::printf("FAIL MC within 3 SE: |%.4f-%.4f|=%.4f > 3*%.4f\n",
                    mc.price, bs, std::fabs(mc.price - bs), mc.std_error);
        ++failures;
    } else {
        std::printf("ok   MC within 3 SE of BS (price %.4f, se %.4f)\n",
                    mc.price, mc.std_error);
    }
    if (mc.std_error <= 0) { std::printf("FAIL MC std_error <= 0\n"); ++failures; }
}

static void test_invalid_option_inputs() {
    OptionSpec base{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};
    auto check_bs = [](const OptionSpec& o, const char* name) {
        check_invalid_argument([&] { (void)bs_price(o); }, name);
    };

    OptionSpec bad = base;
    bad.S = 0.0;
    check_bs(bad, "reject zero spot");
    bad = base; bad.S = -1.0;
    check_bs(bad, "reject negative spot");

    bad = base; bad.K = 0.0;
    check_bs(bad, "reject zero strike");
    bad = base; bad.K = -1.0;
    check_bs(bad, "reject negative strike");

    bad = base; bad.sigma = -0.01;
    check_bs(bad, "reject negative volatility");
    bad = base; bad.sigma = 0.0;
    check_bs(bad, "BS rejects zero volatility");

    bad = base; bad.T = 0.0;
    check_bs(bad, "BS rejects zero maturity");
    bad = base; bad.T = -1.0;
    check_bs(bad, "reject negative maturity");

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    bad = base; bad.S = nan;
    check_bs(bad, "reject NaN spot");
    bad = base; bad.K = inf;
    check_bs(bad, "reject infinite strike");
    bad = base; bad.r = nan;
    check_bs(bad, "reject NaN rate");
    bad = base; bad.q = inf;
    check_bs(bad, "reject infinite dividend yield");
    bad = base; bad.sigma = inf;
    check_bs(bad, "reject infinite volatility");
    bad = base; bad.T = inf;
    check_bs(bad, "reject infinite maturity");

    bad = base; bad.type = static_cast<OptionType>(99);
    check_bs(bad, "reject invalid option type");
    bad = base; bad.S = 0.0;
    check_invalid_argument([&] { (void)bs_delta(bad); },
                           "Greeks reject invalid specifications");
}

static void test_invalid_binomial_inputs() {
    OptionSpec base{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};

    check_invalid_argument([&] { (void)crr_price(base, 0, false); },
                           "CRR rejects zero steps");
    check_invalid_argument([&] { (void)crr_price(base, -1, false); },
                           "CRR rejects negative steps");

    OptionSpec bad = base;
    bad.T = 0.0;
    check_invalid_argument([&] { (void)crr_price(bad, 10, false); },
                           "CRR rejects zero maturity");
    bad = base; bad.sigma = 0.0;
    check_invalid_argument([&] { (void)crr_price(bad, 10, false); },
                           "CRR rejects zero volatility");

    bad = base; bad.r = 1.0; bad.sigma = 0.01;
    check_invalid_argument([&] { (void)crr_price(bad, 1, false); },
                           "CRR rejects probability above one");
    bad = base; bad.q = 1.0; bad.sigma = 0.01;
    check_invalid_argument([&] { (void)crr_price(bad, 1, false); },
                           "CRR rejects probability below zero");
}

static void test_invalid_monte_carlo_inputs() {
    OptionSpec base{100, 100, 0.05, 0.0, 0.20, 1.0, OptionType::Call};

    check_invalid_argument([&] { (void)mc_price(base, 1, 7); },
                           "MC rejects one antithetic pair");
    check_invalid_argument([&] { (void)mc_price(base, 0, 7); },
                           "MC rejects zero antithetic pairs");
    check_invalid_argument([&] { (void)mc_price(base, -1, 7); },
                           "MC rejects negative antithetic pairs");

    OptionSpec bad = base;
    bad.T = -1.0;
    check_invalid_argument([&] { (void)mc_price(bad, 2, 7); },
                           "MC rejects negative maturity");
    bad = base; bad.sigma = -0.01;
    check_invalid_argument([&] { (void)mc_price(bad, 2, 7); },
                           "MC rejects negative volatility");
    bad = base; bad.q = std::numeric_limits<double>::quiet_NaN();
    check_invalid_argument([&] { (void)mc_price(bad, 2, 7); },
                           "MC rejects non-finite inputs");

    OptionSpec expiry = base;
    expiry.S = 110.0;
    expiry.T = 0.0;
    McResult at_expiry = mc_price(expiry, 2, 7);
    CHECK_NEAR(at_expiry.price, 10.0, 1e-12, "MC handles zero maturity");
    CHECK_NEAR(at_expiry.std_error, 0.0, 1e-12, "MC zero-maturity SE");

    OptionSpec deterministic = base;
    deterministic.sigma = 0.0;
    McResult zero_vol = mc_price(deterministic, 2, 7);
    double terminal = deterministic.S
                    * std::exp((deterministic.r - deterministic.q) * deterministic.T);
    double expected = std::exp(-deterministic.r * deterministic.T)
                    * payoff(deterministic.type, terminal, deterministic.K);
    CHECK_NEAR(zero_vol.price, expected, 1e-12, "MC handles zero volatility");
    CHECK_NEAR(zero_vol.std_error, 0.0, 1e-12, "MC zero-volatility SE");
}

int main() {
    test_payoff();
    test_black_scholes();
    test_greeks();
    test_binomial();
    test_monte_carlo();
    test_invalid_option_inputs();
    test_invalid_binomial_inputs();
    test_invalid_monte_carlo_inputs();
    if (failures) { std::printf("\n%d FAILURES\n", failures); return 1; }
    std::printf("\nALL PASS\n");
    return 0;
}
