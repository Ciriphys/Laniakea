#pragma once

#include "lnkpch.h"

#include <Utility/Macro.h>

enum class ComplexForm : int {
	Cartesian = 0,
	Polar = 1
};

struct LNK_API Complex {
	union {
		struct {
			double Re;
			double Im;
		};
		struct {
			double x;
			double y;
		};
		struct {
			double r;
			double w;
		};
	};

	ComplexForm form;

	Complex(double first, double second = 0.0, ComplexForm form = ComplexForm::Cartesian);

	// Operators
	Complex Conjugate();
	double Modulus();
	double Modulus2();
	double Arg();

	// Switch Form
	void ToPolar();
	void ToCartesian();

	// Complex Algebra 
	Complex operator+(Complex other);
	Complex operator-(Complex other);
	Complex operator*(Complex other);
	Complex operator/(Complex other);

	void operator+=(Complex other);
	void operator-=(Complex other);
	void operator*=(Complex other);
	void operator/=(Complex other);

	Complex operator*(double scalar);
	void operator*=(double scalar);

	// Utility
	
	LNK_API friend std::ostream& operator<<(std::ostream& os, Complex& z);
	bool operator==(const Complex& other);
	bool operator!=(const Complex& other);
	Complex operator-();
	Complex operator*();

	static Complex Zero(ComplexForm form = ComplexForm::Cartesian);
	static Complex One(ComplexForm form = ComplexForm::Cartesian);
	static Complex Imaginary(ComplexForm form = ComplexForm::Cartesian);
};

LNK_API void Conjugate(Complex& z);