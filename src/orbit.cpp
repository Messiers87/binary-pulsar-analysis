#include "orbit.h"

// Define physical constants
const long double G = 6.67430e-11L;      // m^3 kg^-1 s^-2
const long double MSUN = 1.98847e30L;    // kg
const long double DAY = 86400.0L;        // seconds
const long double PI = 3.141592653589793238462643383279502884L;
const long double C = 299792458.0L;      // m/s

long double deg2rad(long double d) { return d * PI / 180.0L; }
long double rad2deg(long double r) { return r * 180.0L / PI; }

long double norm2pi(long double x) {
    long double y = fmod(x, 2.0L * PI);
    if (y < 0) y += 2.0L * PI;
    return y;
}

long double solve_kepler_E(long double M, long double e) {
    M = norm2pi(M);
    if (e == 0.0L) return M;

    long double E = (e < 0.8L) ? M : PI;
    for (int it = 0; it < 60; ++it) {
        long double f_val = E - e * sinl(E) - M;
        long double fp = 1.0L - e * cosl(E);
        long double dE = -f_val / fp;
        E += dE;
        if (fabsl(dE) < 1e-14L) break;
    }
    return norm2pi(E);
}

long double E_to_true_f(long double E, long double e) {
    long double s = sqrtl(1.0L + e) * sinl(E / 2.0L);
    long double c = sqrtl(1.0L - e) * cosl(E / 2.0L);
    long double f = 2.0L * atan2l(s, c);
    return norm2pi(f);
}

long double projected_semi_major_axis(const OrbitParams& p) {
    long double Mp = p.Mp_sun * MSUN;
    long double Mc = p.Mc_sun * MSUN;
    long double Po = p.Po_day * DAY;
    long double sin_i = sinl(deg2rad(p.inc_deg));
    long double ap_prime = cbrtl((powl(Po / (2.0L * PI), 2.0L) * G * (Mp + Mc))) * (Mc / (Mp + Mc)) * sin_i;
    return ap_prime;
}

long double r_l(const OrbitParams& p, long double t) {
    long double Po = p.Po_day * DAY;
    long double omega_o = 2.0L * PI / Po;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);
    long double M = omega_o * (t - p.Tp_s);
    long double E = solve_kepler_E(M, e);
    long double f = E_to_true_f(E, e);
    long double varpi = deg2rad(p.varpi_deg);
    long double rl = ap * sinl(f + varpi) * (1.0L - e * e) / (1.0L + e * cosl(f));
    return rl;
}