#pragma once
#include "orbit.h"
#include "precision.h"

real gamma1(const OrbitParams& baseParams, real f, real T_obs, int m, real Pp_seconds);
real compute_w(real f0_deg, real e);

