#pragma once
#include <cmath>

// switch between precision modes
using real = double; // fast mode
// using real = long double; // slow but high-precision mode

// Define mathematical constants
constexpr real PI = static_cast<real>(3.14159265358979323846);
constexpr real C  = static_cast<real>(299792458.0);
constexpr real G  = static_cast<real>(6.67430e-11);
constexpr real MSUN = static_cast<real>(1.98847e30);
constexpr real DAY  = static_cast<real>(86400.0);

// helper overloads for math functions
inline real mycos(real x) { return std::cos(x); }
inline real mysin(real x) { return std::sin(x); }
inline real mypow(real x, real y) { return std::pow(x, y); }
inline real mydeg2rad(real deg) { return deg * static_cast<real>(PI / 180.0); }