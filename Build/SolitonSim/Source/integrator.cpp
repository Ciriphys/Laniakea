#include <cmath>

#include "integrator.hpp"

vec4 soliton_system(vec4 old, double r) 
{
    return  {
            old.df,
            - 2.0 / r * old.df + 2.0 * old.phi * old.f,
            old.dphi,
            4.0 * M_PI * SQR(old.f) - 2.0 / r * old.dphi 
    };
}

vec4 runge_kutta_iv(const vec4& state, double r, double h) 
{
    vec4 k1 = soliton_system(state, r);
    vec4 k2 = soliton_system(state + 0.5 * h * k1, r);
    vec4 k3 = soliton_system(state + 0.5 * h * k2, r);
    vec4 k4 = soliton_system(state + h * k3, r);

    return state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

