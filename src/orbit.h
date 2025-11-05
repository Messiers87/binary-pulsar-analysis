#ifndef ORBIT_H
#define ORBIT_H

#include "precision.h"
#include <cmath>

// ------- Public Declarations of Physical constants -------
extern const real G;
extern const real MSUN;
extern const real DAY;
extern const real PI;
extern const real C;

// ------- Public Declarations of Angle helpers -------
real deg2rad(real d);
real rad2deg(real r);
real norm2pi(real x);

// ------- Public Declarations of Kepler solvers -------
real solve_kepler_E(real M, real e);
real E_to_true_f(real E, real e);

// ------- Public Declaration of the OrbitParams Struct -------
struct OrbitParams {
    real Mp_sun;
    real Mc_sun;
    real Po_day;
    real inc_deg;
    real e;
    real Tp_s;
    real varpi_deg;
};

// ------- Public Declarations of Physics Functions -------
real projected_semi_major_axis(const OrbitParams& p);
real r_l(const OrbitParams& p, real t);
// real v_l(const OrbitParams& p, real t);
// real a_l_of_f(const OrbitParams& p, real f);


#endif // ORBIT_H

