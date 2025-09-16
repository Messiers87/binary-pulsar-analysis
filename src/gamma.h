#pragma once
#include "orbit.h"

long double gamma1(const OrbitParams& baseParams, long double f, long double T_obs, int m, long double Pp_seconds);
long double compute_w(long double f0_deg, long double e);

