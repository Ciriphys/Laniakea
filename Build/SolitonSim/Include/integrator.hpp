#pragma once 

#define SQR(x) ((x) * (x))
 
#include "vec4.hpp" 

vec4 soliton_system(vec4 state, double r);
vec4 runge_kutta_iv(const vec4& state, double r, double h);

