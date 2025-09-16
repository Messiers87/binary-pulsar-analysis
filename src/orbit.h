#pragma once
#include <cmath>

struct OrbitParams {
    long double Mp_sun;
    long double Mc_sun;
    long double Po_day;
    long double inc_deg;
    long double e;
    long double Tp_s;
    long double varpi_deg;
};

long double deg2rad(long double d);
long double rad2deg(long double r);
long double norm2pi(long double x);

long double solve_kepler_E(long double M, long double e);
long double E_to_true_f(long double E, long double e);

long double projected_semi_major_axis(const OrbitParams& p);
long double r_l(const OrbitParams& p, long double t);
long double v_l(const OrbitParams& p, long double t);
long double a_l_of_f(const OrbitParams& p, long double f);
