#include "gamma.h"
#include <vector>
#include <complex>

constexpr long double C = 299792458.0L; 
constexpr long double PI = 3.141592653589793238462643383279502884L;
constexpr long double DAY = 86400.0L; // seconds


// gamma1 function (equation 17) 
long double gamma1(const OrbitParams& baseParams, long double f, long double T_obs, int m, long double Pp_seconds) {

    OrbitParams p = baseParams;

    long double Po = p.Po_day * DAY;
    long double omega_o = 2.0L * PI / Po;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);

    // compute Tp so that f(0) = f0
    long double f0 = deg2rad(f);
    long double tan_halfE = tanl(f0 / 2.0L) * sqrtl((1.0L - e) / (1.0L + e));
    long double E0 = 2.0L * atanl(tan_halfE);

    E0 = norm2pi(E0);
    long double M0 = E0 - e * sinl(E0);
    p.Tp_s = -M0 / omega_o;

    // time sampling
    const int N = 3000; // reduced resolution (was 5000)
    std::vector<long double> t_grid(N + 1);
    std::vector<long double> r_grid(N + 1);
    for (int i = 0; i <= N; i++) {
        t_grid[i] = T_obs * ((long double)i) / ((long double)N);
        r_grid[i] = r_l(p, t_grid[i]);
    }

    long double r0 = r_grid[0];

    // constants
    long double omega_p = 2.0L * PI / Pp_seconds;
    long double K = m * omega_p / C;

    // velocity bounds
    long double vmin = 1.0e10L;
    long double vmax = -1.0e10L;
    for (int i = 1; i <= N; i++) {
        long double v = (r_grid[i] - r_grid[i - 1]) / (t_grid[i] - t_grid[i - 1]);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    long double margin = 1e3L; // add 1 km/s margin

    // helper lambda for trapezoidal integration
    auto compute_gamma = [&](long double alpha) {
        std::complex<long double> sum = 0.0L;
        for (int j = 0; j < N; j++) {
            long double phase1 = K * (r_grid[j] - r0 - alpha * t_grid[j]);
            long double phase2 = K * (r_grid[j + 1] - r0 - alpha * t_grid[j + 1]);
            std::complex<long double> e1 = std::polar(1.0L, phase1);
            std::complex<long double> e2 = std::polar(1.0L, phase2);
            sum += 0.5L * (e1 + e2) * (t_grid[j + 1] - t_grid[j]);
        }
        return std::abs(sum) / T_obs;
        };

    // === Stage 1: coarse scan ===
    int NA_coarse = 100;
    long double alpha_best = 0.0L;
    long double gamma_best = -1.0L;

    for (int i = 0; i < NA_coarse; i++) {
        long double alpha = (vmin - margin) +
            (vmax - vmin + 2.0L * margin) * ((long double)i) / ((long double)NA_coarse - 1);
        long double gamma = compute_gamma(alpha);
        if (gamma > gamma_best) {
            gamma_best = gamma;
            alpha_best = alpha;
        }
    }

    // === Stage 2: refine scan around best alpha ===
    int NA_refine = 100;
    long double refine_width = (vmax - vmin) / 50.0L; // narrower window
    long double alpha_ref_best = alpha_best;
    long double gamma_ref_best = gamma_best;

    for (int i = 0; i < NA_refine; i++) {
        long double alpha = alpha_best - refine_width / 2.0L +
            refine_width * ((long double)i) / ((long double)NA_refine - 1);
        long double gamma = compute_gamma(alpha);
        if (gamma > gamma_ref_best) {
            gamma_ref_best = gamma;
            alpha_ref_best = alpha;
        }
    }

    return gamma_ref_best;
}

// from section 3 
long double compute_w(long double f0_deg, long double e) {
    long double f0 = deg2rad(f0_deg);
    long double num = powl(1.0L + e * cosl(f0), -2.0L);
    long double denom = powl(1.0L + e * cosl(PI), -2.0L);
    return num / denom;
}
