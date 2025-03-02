#include "lnkpch.h"

#include <Linear Algebra/Complex.h>

#include <Utility/Debug.h>

Complex::Complex(double first, double second, ComplexForm form)
{
    this->Re = first;
    this->Im = second;
    this->form = form;
}

Complex Complex::Conjugate()
{
    return {Re, -Im, form};
}

double Complex::Modulus()
{
    return form == ComplexForm::Cartesian ? sqrt(SQR(Re) + SQR(Im)) : r;
}

double Complex::Modulus2()
{
    return form == ComplexForm::Cartesian ? SQR(Re) + SQR(Im) : SQR(r);
}

double Complex::Arg()
{
    return form == ComplexForm::Cartesian ? atan(Im / Re) : w;
}

void Complex::ToPolar()
{
    if (form == ComplexForm::Polar) return;

    double pr = Modulus();
    double pw = Re != 0 ? atan(Im / Re) : Im > 0 ? M_PI / 2 : -M_PI / 2;
    
    form = ComplexForm::Polar;
    r = pr;
    w = pw;
}

void Complex::ToCartesian()
{
    if (form == ComplexForm::Cartesian) return;

    double cx = r * cos(w);
    double cy = r * sin(w);

    form = ComplexForm::Cartesian;
    x = cx;
    y = cy;
}

Complex Complex::operator+(Complex other)
{
    if (form != ComplexForm::Cartesian) ToCartesian();
    if (other.form != ComplexForm::Cartesian) other.ToCartesian();

    return { Re + other.Re, Im + other.Im, ComplexForm::Cartesian };
}

Complex Complex::operator-(Complex other)
{
    if (form != ComplexForm::Cartesian) ToCartesian();
    if (other.form != ComplexForm::Cartesian) other.ToCartesian();

    return { Re - other.Re, Im - other.Im, ComplexForm::Cartesian };
}

Complex Complex::operator*(Complex other)
{
    if (form != ComplexForm::Polar) ToPolar();
    if (other.form != ComplexForm::Polar) other.ToPolar();

    return { r * other.r, w + other.w, ComplexForm::Polar };
}

Complex Complex::operator/(Complex other)
{
    LNK_ASSERT(other != Complex::Zero(), "DomainError: Division by zero!");

    if (form != ComplexForm::Polar) ToPolar();
    if (other.form != ComplexForm::Polar) other.ToPolar();

    return { r / other.r, w - other.w, ComplexForm::Polar };
}

void Complex::operator+=(Complex other)
{
    if (form != ComplexForm::Cartesian) ToCartesian();
    if (other.form != ComplexForm::Cartesian) other.ToCartesian();

    Re += other.Re; Im += other.Im;
}

void Complex::operator-=(Complex other)
{
    if (form != ComplexForm::Cartesian) ToCartesian();
    if (other.form != ComplexForm::Cartesian) other.ToCartesian();

    Re -= other.Re; Im -= other.Im;
}

void Complex::operator*=(Complex other)
{
    if (form != ComplexForm::Polar) ToPolar();
    if (other.form != ComplexForm::Polar) other.ToPolar();

    r *= other.r; w += other.w;
}

void Complex::operator/=(Complex other)
{
    LNK_ASSERT(other != Complex::Zero(), "DomainError: Division by zero!");

    if (form != ComplexForm::Polar) ToPolar();
    if (other.form != ComplexForm::Polar) other.ToPolar();

    r /= other.r; w -= other.w;
}

Complex Complex::operator*(double scalar)
{
    return form == ComplexForm::Cartesian ? Complex(Re * scalar, Im * scalar, form) : Complex(abs(scalar) * r, w, form);
}

void Complex::operator*=(double scalar)
{
    if (form == ComplexForm::Cartesian) {
        Re *= scalar;
        Im *= scalar;
    }
    else {
        r *= abs(scalar);
    }
}

LNK_API std::ostream& operator<<(std::ostream& os, Complex& z)
{
    char sign = z.Im > 0 ? '+' : '-';

    if (z.form == ComplexForm::Cartesian)
    {
        if (z.Im != 0) {
            os << z.Re << ' ' << sign << ' ' << abs(z.Im) << 'i';
        }
        else {
            os << z.Re;
        }
    }
    else
    {
        if (z.r != 0.0) {
            os << z.r << "*exp{" << z.w << "i}";
        }
        else {
            os << z.r;
        }
    }

    return os;
}

bool Complex::operator==(const Complex& other)
{
    return Re == other.Re && Im == other.Im;
}

bool Complex::operator!=(const Complex& other)
{
    return !(operator==(other));
}

Complex Complex::operator-()
{
    if (form != ComplexForm::Cartesian) ToCartesian();
    return {-Re, -Im, ComplexForm::Cartesian};
}

Complex Complex::operator*()
{
    return Conjugate();
}

Complex Complex::Zero(ComplexForm form)
{
    return {0.0, 0.0, form};
}

Complex Complex::One(ComplexForm form)
{
    return {1.0, 0.0, form};
}

Complex Complex::Imaginary(ComplexForm form)
{
    return form == ComplexForm::Cartesian ? Complex(0.0, 1.0, form) : Complex(1.0, M_PI / 2, form);
}

LNK_API void Conjugate(Complex& z)
{
    z = *z;
}
