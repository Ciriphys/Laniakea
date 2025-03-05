#include "lnkpch.h"

#include <Linear Algebra/Complex.h>
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

bool Vec2::operator==(const Vec2& other) 
{
	return x == other.x && y == other.y;
}

bool Vec2::operator!=(const Vec2& other) 
{
	return !operator==(other);
}

Complex Vec2::ToComplex() 
{
	return Complex(x, y);
}

Vec2 Vec2::Zero()
{
	return Vec2(0.0, 0.0);
}

Vec2 Vec2::e1()
{
	return Vec2(1.0, 0.0);
}

Vec2 Vec2::e2()
{
	return Vec2(0.0, 1.0);
}

ComplexVec2::ComplexVec2(Complex _x, Complex _y) : x(_x), y(_y) {}

double ComplexVec2::Norm()
{
	return sqrt(x.Modulus2() + y.Modulus2());
}

double ComplexVec2::Norm2()
{
	return x.Modulus2() + y.Modulus2();
}

ComplexVec2 ComplexVec2::operator+(const ComplexVec2& other)
{
	return ComplexVec2(x + other.x, y + other.y);
}

ComplexVec2 ComplexVec2::operator-(const ComplexVec2& other)
{
	return ComplexVec2(x - other.x, y - other.y);;
}

ComplexVec2 ComplexVec2::operator*(Complex lambda)
{
	return ComplexVec2(x * lambda, y * lambda);
}

// Computes the standard scalar product
Complex ComplexVec2::operator*(ComplexVec2 other)
{
	return x * *other.x + y * *other.y;
}

void ComplexVec2::operator+=(const ComplexVec2& other)
{
	x += other.x;
	y += other.y;
}

void ComplexVec2::operator-=(const ComplexVec2& other)
{
	x -= other.x;
	y -= other.y;
}

void ComplexVec2::operator*=(Complex lambda)
{
	x *= lambda;
	y *= lambda;
}

LNK_API std::ostream& operator<<(std::ostream& os, const ComplexVec2& vector)
{
	os << "(" << vector.x << ", " << vector.y << ")";
	return os;
}

bool ComplexVec2::operator==(const ComplexVec2& other)
{
	return x == other.x && y == other.y;
}

bool ComplexVec2::operator!=(const ComplexVec2& other)
{
	return !operator==(other);
}

ComplexVec2 ComplexVec2::Zero()
{
	return ComplexVec2(0.0, 0.0);
}

ComplexVec2 ComplexVec2::e1()
{
	return ComplexVec2(1.0, 0.0);
}

ComplexVec2 ComplexVec2::e2()
{
	return ComplexVec2(0.0, 1.0);
}