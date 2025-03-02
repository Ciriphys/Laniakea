#include "Laniakea.h"

#include <iostream>

int main() 
{
	Laniakea::DisplayLog("Hello, World!\n");
	Complex z = Complex(2.0, 3.0);
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
	std::cout << z / zero << std::endl;

	return 0;
}