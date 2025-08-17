#pragma once

#include <iostream>

struct vec3 {
    vec3(double x = 0, double y = 0, double z = 0);

    vec3 operator+(const vec3& other) const;
    vec3 operator-(const vec3& other) const;
    vec3 operator*(double scalar) const;
    vec3 operator/(double scalar) const;
    
    double operator*(const vec3& other) const;

    vec3& operator+=(const vec3& other);
    vec3& operator-=(const vec3& other);
    vec3& operator*=(double scalar);
    vec3& operator/=(double scalar);

    friend vec3 operator*(double scalar, const vec3& v);
    friend vec3 operator/(double scalar, const vec3& v);

    friend std::ostream& operator<<(std::ostream& os, const vec3& v);

     union {
        double data[3];
        struct {
            double x, y, z;
        };
    };
}; 