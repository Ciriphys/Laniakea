#pragma once

#include "lnkpch.h"

#include <Linear Algebra/Complex.h>
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

	// Utility
	LNK_API friend std::ostream& operator<<(std::ostream& os, const Vec2& vector);
	bool operator==(const Vec2& other);
	bool operator!=(const Vec2& other);
	Complex ToComplex();

	static Vec2 Zero();
	static Vec2 e1();
	static Vec2 e2();
};

struct LNK_API ComplexVec2 {
	Complex x;
	Complex y;

	ComplexVec2(Complex _x, Complex _y);

	// Functions
	double Norm();
	double Norm2();

	// Vector Algebra 
	ComplexVec2 operator+(const ComplexVec2& other);
	ComplexVec2 operator-(const ComplexVec2& other);
	ComplexVec2 operator*(Complex lambda);
	Complex operator*(ComplexVec2 other);

	void operator+=(const ComplexVec2& other);
	void operator-=(const ComplexVec2& other);
	void operator*=(Complex lambda);

	// Utility
	LNK_API friend std::ostream& operator<<(std::ostream& os, const ComplexVec2& vector);
	bool operator==(const ComplexVec2& other);
	bool operator!=(const ComplexVec2& other);

	static ComplexVec2 Zero();
	static ComplexVec2 e1();
	static ComplexVec2 e2();
};