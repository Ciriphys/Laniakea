#include "lnkpch.h"

#include <Linear Algebra/Vec2.h>
#include <Utility/Debug.h>

Vec2::Vec2(double _x, double _y) : x(_x), y(_y) {}

double Vec2::Norm()
{
	return sqrt(SQR(x) + SQR(y));
}

double Vec2::Norm2()
{
	return SQR(x) + SQR(y);
}

Vec2 Vec2::operator+(const Vec2& other)
{
	return Vec2(x + other.x, y + other.y);
}

Vec2 Vec2::operator-(const Vec2& other)
{
	return Vec2(x - other.x, y - other.y);;
}

Vec2 Vec2::operator*(double lambda)
{
	return Vec2(lambda * x, lambda * y);
}

// Computes the standard scalar product
double Vec2::operator*(const Vec2& other)
{
	return x * other.x + y * other.y;
}

void Vec2::operator+=(const Vec2& other)
{
	x += other.x;
	y += other.y;
}

void Vec2::operator-=(const Vec2& other)
{
	x -= other.x;
	y -= other.y;
}

void Vec2::operator*=(double lambda)
{
	x *= lambda;
	y *= lambda;
}

LNK_API std::ostream& operator<<(std::ostream& os, const Vec2& vector)
{
	os << "(" << vector.x << ", " << vector.y << ")";
	return os;
}