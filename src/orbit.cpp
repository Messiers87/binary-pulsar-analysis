#include "orbit.h"
#include "precision.h"

real deg2rad(real d) { return d * PI / static_cast<real>(180.0); }
real rad2deg(real r) { return r * static_cast<real>(180.0) / PI; }

real norm2pi(real x) {
    real y = fmod(x, static_cast<real>(2.0) * PI);
    if (y < 0) y += static_cast<real>(2.0) * PI;
    return y;
}

real solve_kepler_E(real M, real e) {
    M = norm2pi(M);
    if (e == static_cast<real>(0.0)) return M;

    real E = (e < static_cast<real>(0.8)) ? M : PI;
    for (int it = 0; it < 60; ++it) {
        real f_val = E - e * sinl(E) - M;
        real fp = static_cast<real>(1.0) - e * cosl(E);
        real dE = -f_val / fp;
        E += dE;
        if (fabsl(dE) < 1e-14L) break;
    }
    return norm2pi(E);
}

real E_to_true_f(real E, real e) {
    real s = std::sqrt(static_cast<real>(1.0) + e) * sinl(E / static_cast<real>(2.0));
    real c = std::sqrt(static_cast<real>(1.0) - e) * cosl(E / static_cast<real>(2.0));
    real f = static_cast<real>(2.0) * atan2l(s, c);
    return norm2pi(f);
}

real projected_semi_major_axis(const OrbitParams& p) {
    real Mp = p.Mp_sun * MSUN;
    real Mc = p.Mc_sun * MSUN;
    real Po = p.Po_day * DAY;
    real aR = std::cbrt((std::pow(Po / (static_cast<real>(2.0) * PI), static_cast<real>(2.0)) * G * (Mp + Mc)));
    real sin_i = std::sin(deg2rad(p.inc_deg));

    real ap_prime = aR * (Mc / (Mp + Mc)) * sin_i; // meters
    return ap_prime;
}

real r_l(const OrbitParams& p, real t) {
    real Po = p.Po_day * DAY;
    real omega_o = static_cast<real>(2.0) * PI / Po;
    real e = p.e;
    real ap = projected_semi_major_axis(p);
    real M = omega_o * (t - p.Tp_s);
    real E = solve_kepler_E(M, e);
    real f = E_to_true_f(E, e);
    real varpi = deg2rad(p.varpi_deg);

    real rl = ap * sinl(f + varpi) * (static_cast<real>(1.0) - e * e) / (static_cast<real>(1.0) + e * cosl(f));
    return rl;
}