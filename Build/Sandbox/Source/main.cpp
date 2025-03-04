#include "Laniakea.h"

#include <iostream>

void TestComplex()
{
	Complex z = Complex(2.0, 3.0, ComplexForm::Cartesian);
	Complex w = Complex::One();

	std::cout << z + w << std::endl;
	std::cout << z * w << std::endl;

	z += w;
	std::cout << z.Conjugate() << std::endl;
	std::cout << z + *z << std::endl;

	z.ToPolar();
	std::cout << z << std::endl;

	Complex zero = Complex::Zero();
	std::cout << zero * 5.0 << std::endl;
	std::cout << z * zero << std::endl;
	std::cout << z + zero << std::endl;
	// Assertion! std::cout << z / zero << std::endl;
}

void TestVec2() 
{
	Vec2 u = {  3.0, 1.0 };
	Vec2 v = { -1.0, 3.0 };

	std::cout << u.Norm() << " " << v.Norm2() << std::endl;
	std::cout << u << " " << v << " " << u + v << std::endl;
	std::cout << (u + v) * (u * v) << std::endl;
	std::cout << u - v << std::endl;
	u -= v;
	std::cout << u << std::endl;
}

int main() 
{
	Laniakea::DisplayLog("Hello, World!\n");
	TestVec2();

	return 0;
}