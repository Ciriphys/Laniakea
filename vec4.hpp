#pragma once

#include <iostream>

struct vec4 {
    vec4(double x = 0, double y = 0, double z = 0, double t = 0);

    vec4 operator+(const vec4& other) const;
    vec4 operator-(const vec4& other) const;

    vec4 operator*(double scalar) const;
    vec4 operator/(double scalar) const;

    vec4& operator+=(const vec4& other);
    vec4& operator-=(const vec4& other);
    vec4& operator*=(double scalar);
    vec4& operator/=(double scalar);

    friend vec4 operator*(double scalar, const vec4& v);
    friend vec4 operator/(double scalar, const vec4& v);

    friend std::ostream& operator<<(std::ostream& os, const vec4& v);

    union {
        struct {
            double x, y, z, t;
        };
        struct {
            double f, df, phi, dphi;
        };
    }; 
}; 