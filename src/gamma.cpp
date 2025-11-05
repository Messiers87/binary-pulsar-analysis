#include "gamma.h"
#include "../src/precision.h"
#include <vector>
#include <complex>

// gamma1 function (equation 17) 
real gamma1(const OrbitParams& baseParams, real f, real T_obs, int m, real Pp_seconds) {

    OrbitParams p = baseParams;

    real Po = p.Po_day * DAY;
    real omega_o = static_cast<real>(2.0)* PI / Po;
    real e = p.e;
    real ap = projected_semi_major_axis(p);

    // compute Tp so that f(0) = f0
    real f0 = deg2rad(f);
    real tan_halfE = std::tan(f0 / static_cast<real>(2.0)) * std::sqrt((static_cast<real>(1.0) - e) / (static_cast<real>(1.0) + e));
    real E0 = static_cast<real>(2.0) * std::atan(tan_halfE);

    E0 = norm2pi(E0);
    real M0 = E0 - e * sinl(E0);
    p.Tp_s = -M0 / omega_o;

    // time sampling
    const int N = 3000; // resolution
    std::vector<real> t_grid(N + 1);
    std::vector<real> r_grid(N + 1);
    for (int i = 0; i <= N; i++) {
        t_grid[i] = T_obs * ((real)i) / ((real)N);
        r_grid[i] = r_l(p, t_grid[i]);
    }

    real r0 = r_grid[0];

    // constants
    real omega_p = static_cast<real>(2.0) * PI / Pp_seconds;
    real K = static_cast<real>(m) * omega_p / C;

    // velocity bounds
    real vmin = static_cast<real>(1.0e10);
    real vmax = static_cast<real>(-1.0e10);
    for (int i = 1; i <= N; i++) {
        real v = (r_grid[i] - r_grid[i - 1]) / (t_grid[i] - t_grid[i - 1]);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }

    real margin = std::fabs(vmax - vmin) * static_cast<real>(1.0) + static_cast<real>(1e2);

    // helper lambda for trapezoidal integration
    auto compute_gamma = [&](real alpha) {
        std::complex<real> sum = static_cast<real>(0.0);
        for (int j = 0; j < N; j++) {
            real phase1 = K * (r_grid[j] - r0 - alpha * t_grid[j]);
            real phase2 = K * (r_grid[j + 1] - r0 - alpha * t_grid[j + 1]);
            std::complex<real> e1 = std::polar(static_cast<real>(1.0), phase1);
            std::complex<real> e2 = std::polar(static_cast<real>(1.0), phase2);
            sum += static_cast<real>(0.5) * (e1 + e2) * (t_grid[j + 1] - t_grid[j]);
        }
        return std::abs(sum) / T_obs;
        };

    // === Stage 1: coarse scan ===
    int NA_coarse = 100;
    real alpha_best = static_cast<real>(0.0);
    real gamma_best = static_cast<real>(-1.0);

    for (int i = 0; i < NA_coarse; i++) {
        real alpha = (vmin - margin) +
            (vmax - vmin + static_cast<real>(2.0) * margin) * (static_cast<real>(i)) / (static_cast<real>(NA_coarse) - static_cast<real>(1.0));
        real gamma = compute_gamma(alpha);
        if (gamma > gamma_best) {
            gamma_best = gamma;
            alpha_best = alpha;
        }
    }

    // === Stage 2: refine scan around best alpha ===
    int NA_refine = 200;
	real coarse_stepsize = (vmax - vmin + static_cast<real>(2.0) * margin) / ((static_cast<real>(NA_coarse) - static_cast<real>(1.0)));
	real refine_width = static_cast<real>(2.0) * coarse_stepsize; // twice the coarse stepsize
    // real refine_width = (vmax - vmin) / 50.0L; // narrower window

    real alpha_ref_best = alpha_best;
    real gamma_ref_best = gamma_best;

    for (int i = 0; i < NA_refine; i++) {
        real alpha = alpha_best - refine_width / static_cast<real>(2.0) +
            refine_width * (static_cast<real>(i)) / (static_cast<real>(NA_refine) - static_cast<real>(1.0));
        real gamma = compute_gamma(alpha);
        if (gamma > gamma_ref_best) {
            gamma_ref_best = gamma;
            alpha_ref_best = alpha;
        }
    }

    return gamma_ref_best;
}

// from section 3 
real compute_w(real f0_deg, real e) {
    real f0 = mydeg2rad(f0_deg);
    real num = mypow(static_cast<real>(1.0) + e * cosl(f0), static_cast<real>(-2.0));
    real denom = mypow(static_cast<real>(1.0) + e * cosl(PI), static_cast<real>(-2.0));
    return num / denom;
}
