#include "vec3.hpp"

vec3::vec3(double x, double y, double z)
    : x(x), y(y), z(z) {}

vec3 vec3::operator+(const vec3& other) const {
    return vec3(x + other.x, y + other.y, z + other.z);
}

vec3 vec3::operator-(const vec3& other) const {
    return vec3(x - other.x, y - other.y, z - other.z);
}

vec3 vec3::operator*(double scalar) const {
    return vec3(x * scalar, y * scalar, z * scalar);
}

vec3 vec3::operator/(double scalar) const {
    return vec3(x / scalar, y / scalar, z / scalar);
}

vec3& vec3::operator+=(const vec3& other) {
    x += other.x; y += other.y; z += other.z;
    return *this;
}

vec3& vec3::operator-=(const vec3& other) {
    x -= other.x; y -= other.y; z -= other.z;
    return *this;
}

vec3& vec3::operator*=(double scalar) {
    x *= scalar; y *= scalar; z *= scalar;
    return *this;
}

vec3& vec3::operator/=(double scalar) {
    x /= scalar; y /= scalar; z /= scalar;
    return *this;
}

vec3 operator*(double scalar, const vec3& v) {
    return v * scalar;
}

double vec3::operator*(const vec3& other) const {
    return x * other.x + y * other.y + z * other.z;
}

vec3 operator/(double scalar, const vec3& v) {
    return vec3(scalar / v.x, scalar / v.y, scalar / v.z);
}

std::ostream& operator<<(std::ostream& os, const vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}