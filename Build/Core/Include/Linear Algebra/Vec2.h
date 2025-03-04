#pragma once

#include "lnkpch.h"

#include <Utility/Macro.h>

struct LNK_API Vec2 {
	double x;
	double y;

	Vec2(double _x, double _y);
	
	// Functions
	double Norm();
	double Norm2();

	// Vector Algebra 
	Vec2 operator+(const Vec2& other);
	Vec2 operator-(const Vec2& other);
	Vec2 operator*(double lambda);
	double operator*(const Vec2& other);

	void operator+=(const Vec2& other);
	void operator-=(const Vec2& other);
	void operator*=(double lambda);

	LNK_API friend std::ostream& operator<<(std::ostream& os, const Vec2& vector);
};