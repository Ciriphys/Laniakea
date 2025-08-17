#include "vec4.hpp"

vec4::vec4(double x, double y, double z, double t)
    : x(x), y(y), z(z), t(t) {}

vec4 vec4::operator+(const vec4& other) const {
    return vec4(x + other.x, y + other.y, z + other.z, t + other.t);
}

vec4 vec4::operator-(const vec4& other) const {
    return vec4(x - other.x, y - other.y, z - other.z, t - other.t);
}

vec4 vec4::operator*(double scalar) const {
    return vec4(x * scalar, y * scalar, z * scalar, t * scalar);
}

vec4 vec4::operator/(double scalar) const {
    return vec4(x / scalar, y / scalar, z / scalar, t / scalar);
}

vec4& vec4::operator+=(const vec4& other) {
    x += other.x; y += other.y; z += other.z; t += other.t;
    return *this;
}

vec4& vec4::operator-=(const vec4& other) {
    x -= other.x; y -= other.y; z -= other.z; t -= other.t;
    return *this;
}

vec4& vec4::operator*=(double scalar) {
    x *= scalar; y *= scalar; z *= scalar; t *= scalar;
    return *this;
}

vec4& vec4::operator/=(double scalar) {
    x /= scalar; y /= scalar; z /= scalar; t /= scalar;
    return *this;
}

vec4 operator*(double scalar, const vec4& v) {
    return v * scalar;
}

vec4 operator/(double scalar, const vec4& v) {
    return vec4(scalar / v.x, scalar / v.y, scalar / v.z, scalar / v.t);
}

std::ostream& operator<<(std::ostream& os, const vec4& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.t << ")";
    return os;
}