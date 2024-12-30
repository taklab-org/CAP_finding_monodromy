/*
 * variable_coefficients.hpp:
 * Define the vector field of the Maurer-Cartan form
 *
 * written by A. Takayasu	Jan. 30, 2020
 */

#ifndef VARIABLE_COEFFICIENTS_HPP
#define VARIABLE_COEFFICIENTS_HPP

// #include <kv/complex.hpp>
// #include <kv/interval.hpp>
// #include <kv/dd.hpp>

// typedef double dd;
// typedef kv::dd dd;
// typedef kv::interval<dd> intval;
namespace ub = boost::numeric::ublas;

namespace kv {
	// h(x,y)
	template <class T> T h(const T x, const T y)
	{
		return (T) 1. + (T) 20. * x + (T) 9. * y;
	}
	// a(x,y)
	template <class T> T a(const T x, const T y)
	{
		T A1 = (T) 4. * x + (T) 16. * pow(x,2) - (T) 3. * y - (T) 60. * x * y - (T) 27. * pow(y,2);
		return A1 / ((T) 2. * x * y * h(x,y));
	}
	// a_y(x,y)
	template <class T> T a_y(const T x, const T y)
	{
		// T hy = (T) 9.;
		// T A1 = ((T) 4. * x + (T) 16. * pow(x,2) + (T) 27. * pow(y,2)) / ((T) 2. * x * pow(y,2) * h(x,y));
		// T A2 = ((T) 4. * x + (T) 16. * pow(x,2) - (T) 3. * y - (T) 60. * x * y - (T) 27. * pow(y,2)) / ((T) 2. * x * y * pow(h(x,y),2)) * hy;
		// return -A1 - A2;
		return -((T) 2. * ((T) 4. * x + (T) 1.) * ((T) 20. * x + (T) 18. * y + (T) 1.)) / (pow(y,2) * pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// b(x,y)
	template <class T> T b(const T x, const T y)
	{
		T A1 = (T) 16. * x + (T) 96. * pow(x,2) + (T) 4. * y + (T) 168. * x * y + (T) 27. * pow(y,2);
		return -A1 / ((T) 4. * pow(x,2) * h(x,y));
	}
	// b_y(x,y)
	template <class T> T b_y(const T x, const T y)
	{
		T hy = (T) 9.;
		T A1 = ((T) 4. + (T) 168. * x + (T) 54. * y) / ((T) 4. * pow(x,2) * h(x,y));
		T A2 = ((T) 16. * x + (T) 96. * pow(x,2) + (T) 4. * y + (T) 168. * x * y + (T) 27. * pow(y,2)) / ((T) 4. * pow(x,2) * pow(h(x,y),2)) * hy;
		return -A1 + A2;
		// return -((T) 2496. * pow(x,2) + (T) 1080. * x * y + (T) 104. * x + (T) 243. * pow(y,2) + (T) 54. * y + (T) 4.) / ((T) 4. * pow(x,2) * pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// c(x,y)
	template <class T> T c(const T x, const T y)
	{
		return x * ((T) 1. + (T) 4. * x) / (pow(y,2) * h(x,y));
	}
	// c_x(x,y)
	template <class T> T c_x(const T x, const T y)
	{
		// T hx = (T) 20.;
		// T A1 = ((T) 1. + (T) 8. * x) / (pow(y,2) * h(x,y));
		// T A2 = (x * ((T) 1. + (T) 4. * x)) / (pow(y,2) * pow(h(x,y),2)) * hx;
		// return A1 - A2;
		return ((T) 8. * x + (T) 9. * y + (T) 72. * x *y + (T) 80. * pow(x,2) + (T) 1.) / (pow(y,2)*pow((T) 20. * x + (T) 9. * y + 1,2));
	}
	// d(x,y)
	template <class T> T d(const T x, const T y)
	{
		return -((T) 12. * x + (T) 16. * pow(x,2) + y + (T) 72. * x * y) / ((T) 8. * x * y * h(x,y));
	}
	// d_x(x,y)
	template <class T> T d_x(const T x, const T y)
	{
		// T hx = (T) 20.;
		// T A1 = ((T) 16. * pow(x,2) - y) / ((T) 8. * pow(x,2) * y * h(x,y));
		// T A2 = ((T) 12. * x + (T) 16. * pow(x,2) + y + (T) 72. * x * y) / ((T) 8. * x * y * pow(h(x,y),2)) * hx;
		// return -A1 + A2;
		return ((T) 1296. * pow(x,2) * y + (T) 224. * pow(x,2) + (T) 40. * x * y + (T) 9. * pow(y,2) + y) / ((T) 8. * pow(x,2) * y * pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// ell(x,y)
	template <class T> T ell(const T x, const T y)
	{
		T A1 = (T) 8. * x + (T) 32. * pow(x,2) + (T) 4. * y + (T) 84. * x * y + (T) 27. * pow(y,2);
		return -A1 / ((T) 2. * x * h(x,y));
	}
	// ell_y(x,y)
	template <class T> T ell_y(const T x, const T y)
	{
		T hy = (T) 9.;
		T A1 = ((T) 4. + (T) 84. * x + (T) 54. * y) / ((T) 2. * x * h(x,y));
		T A2 = ((T) 8. * x + (T) 32. * pow(x,2) + (T) 4. * y + (T) 84. * x * y + (T) 27. * pow(y,2)) / ((T) 2. * x * pow(h(x,y),2)) * hy;
		return -A1 + A2;
		// return -((T) 1392. * pow(x,2) + (T) 1080. * x * y + (T) 92. * x + (T) 243. * pow(y,2) + (T) 54. * y + (T) 4.) / ((T) 2. * x * pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// em(x,y)
	template <class T> T em(const T x, const T y)
	{
		T A1 = (T) 8. * x + (T) 32. * pow(x,2) + y + (T) 24. * x *y;
		return -A1/((T) 4. * y * h(x,y));
	}
	// em_x(x,y)
	template <class T> T em_x(const T x, const T y)
	{
		T hx = (T) 20.;
		T A1 = ((T) 2. + (T) 16. * x + (T) 6. * y) / ( y * h(x,y));
		T A2 = ((T) 8. * x + (T) 32. * pow(x,2) + y + (T) 24. * x * y) / ((T) 4. * y * pow(h(x,y),2)) * hx;
		return -A1 + A2;
		// return -((T) 160. * pow(x,2) + (T) 144. * x * y + (T) 16. * x + (T) 54. * pow(y,2) + (T) 19. * y + (T) 2.) / (y * pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// p(x,y)
	template <class T> T p(const T x, const T y)
	{
		return ((T) 2. + (T) 12. * x + (T) 9. * y) / (x * y *h(x,y));
	}
	// p_y(x,y)
	template <class T> T p_y(const T x, const T y)
	{
		// T hy = (T) 9.;
		// T A1 = ((T) 2. + (T) 12. * x) / (x * pow(y,2) * h(x,y));
		// T A2 = ((T) 2. + (T) 12. * x + (T) 9. * y) / (x * y * pow(h(x,y),2)) * hy;
		// return -A1 - A2;
		return -(240*pow(x,2) + 216*x*y + 52*x + 81*pow(y,2) + 36*y + 2)/(x*pow(y,2)*pow(h(x,y),2));
	}
	// q(x,y)
	template <class T> T q(const T x, const T y)
	{
		return ((T) 1. - (T) 8. * x) / ((T) 2 * pow(y,2) * h(x,y));
	}
	// q_x(x,y)
	template <class T> T q_x(const T x, const T y)
	{
		T hx = (T) 20.;
		T A1 = (T) 8. / ((T) 2. * pow(y,2) * h(x,y));
		T A2 = ((T) 1. - (T) 8. * x) / ((T) 2 * pow(y,2) * pow(h(x,y),2)) * hx;
		return -A1 - A2;
		// return -((T) 2. * ((T) 18. * y + (T) 7.)) / (pow(y,2)*pow((T) 20. * x + (T) 9. * y + (T) 1.,2));
	}
	// B0(x,y)
	template <class T> T B0(const T x, const T y)
	{
		T A1 = p_y(x,y) + b(x,y) * q(x,y) + ell(x,y) * (q_x(x,y) + c(x,y) * p(x,y));
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return -((((36*y + 14)/(pow(y,2)*pow(h(x,y),2)) - ((4*x + 1)*(12*x + 9*y + 2))/(pow(y,3)*pow(h(x,y),2)))*(32*pow(x,2) + 84*x*y + 8*x + 27*pow(y,2) + 4*y))/(2*x*(20*x + 9*y + 1)) - (240*pow(x,2) + 216*x*y + 52*x + 81*pow(y,2) + 36*y + 2)/(x*pow(y,2)*pow(h(x,y),2)) + ((8*x - 1)*(96*pow(x,2) + 168*x*y + 16*x + 27*pow(y,2) + 4*y))/(8*pow(x,2)*pow(y,2)*pow(h(x,y),2)))/(((8*x + y + 24*x*y + 32*pow(x,2))*(32*pow(x,2) + 84*x*y + 8*x + 27*pow(y,2) + 4*y))/(8*x*y*pow(h(x,y),2)) - 1);
		// return (6144*pow(x,5) + 43776*pow(x,4)*y + 4096*pow(x,4) + 30720*pow(x,3)*pow(y,2) + 16832*pow(x,3)*y + 896*pow(x,3) + 3888*pow(x,2)*pow(y,3) + 12224*pow(x,2)*pow(y,2) + 1856*pow(x,2)*y + 64*pow(x,2) + 3672*x*pow(y,3) + 928*x*pow(y,2) + 64*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2))/(x*pow(y,2)*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// B1(x,y)
	template <class T> T B1(const T x, const T y)
	{
		T A1 = a_y(x,y) + b(x,y) * c(x,y) + ell(x,y) * (c_x(x,y) + c(x,y) * a(x,y)) + ell(x,y) * q(x,y);
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return (2*(2048*pow(x,5) + 18176*pow(x,4)*y + 1536*pow(x,4) + 14592*pow(x,3)*pow(y,2) + 7680*pow(x,3)*y + 384*pow(x,3) + 3888*pow(x,2)*pow(y,3) + 6896*pow(x,2)*pow(y,2) + 944*pow(x,2)*y + 32*pow(x,2) + 2268*x*pow(y,3) + 556*x*pow(y,2) + 40*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)))/(pow(y,2)*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return (2*(2048*pow(x,5) + 18176*pow(x,4)*y + 1536*pow(x,4) + 14592*pow(x,3)*pow(y,2) + 7680*pow(x,3)*y + 384*pow(x,3) + 3888*pow(x,2)*pow(y,3) + 6896*pow(x,2)*pow(y,2) + 944*pow(x,2)*y + 32*pow(x,2) + 2268*x*pow(y,3) + 556*x*pow(y,2) + 40*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)))/(pow(y,2)*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// B2(x,y)
	template <class T> T B2(const T x, const T y)
	{
		T A1 = b_y(x,y) + b(x,y) * d(x,y) + ell(x,y) * (d_x(x,y) + b(x,y) * c(x,y)) + p(x,y);
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return -(49152*pow(x,6) + 399360*pow(x,5)*y + 32768*pow(x,5) + 258048*pow(x,4)*pow(y,2) + 161280*pow(x,4)*y + 7168*pow(x,4) + 10368*pow(x,3)*pow(y,3) + 129536*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 512*pow(x,3) + 44592*pow(x,2)*pow(y,3) + 10688*pow(x,2)*pow(y,2) + 640*pow(x,2)*y + 2916*x*pow(y,4) + 948*x*pow(y,3) + 80*x*pow(y,2) - 243*pow(y,5) - 63*pow(y,4) - 4*pow(y,3))/(4*pow(x,2)*y*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return -(49152*pow(x,6) + 399360*pow(x,5)*y + 32768*pow(x,5) + 258048*pow(x,4)*pow(y,2) + 161280*pow(x,4)*y + 7168*pow(x,4) + 10368*pow(x,3)*pow(y,3) + 129536*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 512*pow(x,3) + 44592*pow(x,2)*pow(y,3) + 10688*pow(x,2)*pow(y,2) + 640*pow(x,2)*y + 2916*x*pow(y,4) + 948*x*pow(y,3) + 80*x*pow(y,2) - 243*pow(y,5) - 63*pow(y,4) - 4*pow(y,3))/(4*pow(x,2)*y*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// B3(x,y)
	template <class T> T B3(const T x, const T y)
	{
		T A1 = ell_y(x,y) + a(x,y) + b(x,y) * em(x,y) + ell(x,y) * (em_x(x,y) + d(x,y) + c(x,y) * ell(x,y));
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return -(8192*pow(x,6) + 124928*pow(x,5)*y + 6144*pow(x,5) + 87552*pow(x,4)*pow(y,2) + 59392*pow(x,4)*y + 1536*pow(x,4) + 12672*pow(x,3)*pow(y,3) + 71008*pow(x,3)*pow(y,2) + 8320*pow(x,3)*y + 128*pow(x,3) + 33120*pow(x,2)*pow(y,3) + 6616*pow(x,2)*pow(y,2) + 320*pow(x,2)*y + 5670*x*pow(y,4) + 1404*x*pow(y,3) + 88*x*pow(y,2) + 243*pow(y,5) + 63*pow(y,4) + 4*pow(y,3))/(x*y*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return -(8192*pow(x,6) + 124928*pow(x,5)*y + 6144*pow(x,5) + 87552*pow(x,4)*pow(y,2) + 59392*pow(x,4)*y + 1536*pow(x,4) + 12672*pow(x,3)*pow(y,3) + 71008*pow(x,3)*pow(y,2) + 8320*pow(x,3)*y + 128*pow(x,3) + 33120*pow(x,2)*pow(y,3) + 6616*pow(x,2)*pow(y,2) + 320*pow(x,2)*y + 5670*x*pow(y,4) + 1404*x*pow(y,3) + 88*x*pow(y,2) + 243*pow(y,5) + 63*pow(y,4) + 4*pow(y,3))/(x*y*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// matA
	template <class T> bool matA(const T x, const T y, ub::matrix<T>& A)
	{
		A(0,0) = (T) 0.0; A(0,1) = (T) 1.0;	A(0,2) = (T) 0.0;	A(0,3) = (T) 0.0;
		A(1,0) = p(x,y);  A(1,1) = a(x,y);  A(1,2) = b(x,y);  A(1,3) = ell(x,y);
		A(2,0) = (T) 0.0; A(2,1) = (T) 0.0;	A(2,2) = (T) 0.0;	A(2,3) = (T) 1.0;
		A(3,0) = B0(x,y); A(3,1) = B1(x,y); A(3,2) = B2(x,y); A(3,3) = B3(x,y);
		// std::cout << "B0: " << B0(x,y) << '\n';
		// std::cout << "B1: " << B1(x,y) << '\n';
		// std::cout << "B2: " << B2(x,y) << '\n';
		// std::cout << "B3: " << B3(x,y) << '\n';
		return true;
	}
	// C0(x,y)
	template <class T> T C0(const T x, const T y)
	{
		T A1 = q_x(x,y) + c(x,y) * p(x,y) + em(x,y) * (p_y(x,y) + b(x,y) * q(x,y));
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return -(67584*pow(x,5) + 76800*pow(x,4)*y + 35840*pow(x,4) + 10368*pow(x,3)*pow(y,2) + 30720*pow(x,3)*y + 6016*pow(x,3) + 12288*pow(x,2)*pow(y,2) + 5088*pow(x,2)*y + 320*pow(x,2) + 1080*x*pow(y,3) + 736*x*pow(y,2) + 64*x*y + 27*pow(y,3) + 4*pow(y,2))/(4*x*pow(y,2)*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return -(67584*pow(x,5) + 76800*pow(x,4)*y + 35840*pow(x,4) + 10368*pow(x,3)*pow(y,2) + 30720*pow(x,3)*y + 6016*pow(x,3) + 12288*pow(x,2)*pow(y,2) + 5088*pow(x,2)*y + 320*pow(x,2) + 1080*x*pow(y,3) + 736*x*pow(y,2) + 64*x*y + 27*pow(y,3) + 4*pow(y,2))/(4*x*pow(y,2)*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// C1(x,y)
	template <class T> T C1(const T x, const T y)
	{
		T A1 = c_x(x,y) + a(x,y) * c(x,y) + em(x,y) * (a_y(x,y) + b(x,y) * c(x,y)) + q(x,y);
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return -(43008*pow(x,5) + 30720*pow(x,4)*y + 25088*pow(x,4) + 10368*pow(x,3)*pow(y,2) + 23296*pow(x,3)*y + 4480*pow(x,3) + 7392*pow(x,2)*pow(y,2) + 3616*pow(x,2)*y + 224*pow(x,2) + 756*x*pow(y,3) + 640*x*pow(y,2) + 56*x*y + 27*pow(y,3) + 4*pow(y,2))/(2*pow(y,2)*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return -(43008*pow(x,5) + 30720*pow(x,4)*y + 25088*pow(x,4) + 10368*pow(x,3)*pow(y,2) + 23296*pow(x,3)*y + 4480*pow(x,3) + 7392*pow(x,2)*pow(y,2) + 3616*pow(x,2)*y + 224*pow(x,2) + 756*x*pow(y,3) + 640*x*pow(y,2) + 56*x*y + 27*pow(y,3) + 4*pow(y,2))/(2*pow(y,2)*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// C2(x,y)
	template <class T> T C2(const T x, const T y)
	{
		T A1 = d_x(x,y) + b(x,y) * c(x,y) + em(x,y) * (b_y(x,y) + b(x,y) * d(x,y)) + em(x,y) * p(x,y);
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return (540672*pow(x,6) + 724992*pow(x,5)*y + 286720*pow(x,5) + 27648*pow(x,4)*pow(y,2) + 294912*pow(x,4)*y + 48128*pow(x,4) + 150400*pow(x,3)*pow(y,2) + 46592*pow(x,3)*y + 2560*pow(x,3) + 14832*pow(x,2)*pow(y,3) + 8096*pow(x,2)*pow(y,2) + 640*pow(x,2)*y - 648*x*pow(y,4) + 372*x*pow(y,3) + 48*x*pow(y,2) + 27*pow(y,4) + 4*pow(y,3))/(16*pow(x,2)*y*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return (540672*pow(x,6) + 724992*pow(x,5)*y + 286720*pow(x,5) + 27648*pow(x,4)*pow(y,2) + 294912*pow(x,4)*y + 48128*pow(x,4) + 150400*pow(x,3)*pow(y,2) + 46592*pow(x,3)*y + 2560*pow(x,3) + 14832*pow(x,2)*pow(y,3) + 8096*pow(x,2)*pow(y,2) + 640*pow(x,2)*y - 648*x*pow(y,4) + 372*x*pow(y,3) + 48*x*pow(y,2) + 27*pow(y,4) + 4*pow(y,3))/(16*pow(x,2)*y*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// C3(x,y)
	template <class T> T C3(const T x, const T y)
	{
		T A1 = em_x(x,y) + d(x,y) + c(x,y) * ell(x,y) + em(x,y) * (ell_y(x,y) + a(x,y) + b(x,y) * em(x,y));
		return A1 / ((T) 1. - ell(x,y) * em(x,y));
		// return (262144*pow(x,6) + 262144*pow(x,5)*y + 155648*pow(x,5) + 36864*pow(x,4)*pow(y,2) + 144384*pow(x,4)*y + 28672*pow(x,4) + 45184*pow(x,3)*pow(y,2) + 25600*pow(x,3)*y + 1536*pow(x,3) - 3456*pow(x,2)*pow(y,3) + 4416*pow(x,2)*pow(y,2) + 448*pow(x,2)*y - 1944*x*pow(y,4) + 72*x*pow(y,3) + 32*x*pow(y,2) + 27*pow(y,4) + 4*pow(y,3))/(8*x*y*(20*x + 9*y + 1)*(1024*pow(x,4) + 256*pow(x,3)*y + 512*pow(x,3) + 704*pow(x,2)*y + 64*pow(x,2) + 252*x*pow(y,2) + 32*x*y + 27*pow(y,3) + 4*pow(y,2)));
		// return (262144*pow(x,6) + 262144*pow(x,5)*y + 155648*pow(x,5) + 36864*pow(x,4)*pow(y,2) + 144384*pow(x,4)*y + 28672*pow(x,4) + 45184*pow(x,3)*pow(y,2) + 25600*pow(x,3)*y + 1536*pow(x,3) - 3456*pow(x,2)*pow(y,3) + 4416*pow(x,2)*pow(y,2) + 448*pow(x,2)*y - 1944*x*pow(y,4) + 72*x*pow(y,3) + 32*x*pow(y,2) + 27*pow(y,4) + 4*pow(y,3))/(8*x*y*(20480*pow(x,5) + 14336*pow(x,4)*y + 11264*pow(x,4) + 2304*pow(x,3)*pow(y,2) + 18944*pow(x,3)*y + 1792*pow(x,3) + 11376*pow(x,2)*pow(y,2) + 1920*pow(x,2)*y + 64*pow(x,2) + 2808*x*pow(y,3) + 620*x*pow(y,2) + 32*x*y + 243*pow(y,4) + 63*pow(y,3) + 4*pow(y,2)));
	}
	// matB
	template <class T> bool matB(const T x, const T y, ub::matrix<T>& B)
	{
		B(0,0) = (T) 0.0; B(0,1) = (T) 0.0;	B(0,2) = (T) 1.0;	B(0,3) = (T) 0.0;
		B(1,0) = (T) 0.0; B(1,1) = (T) 0.0;	B(1,2) = (T) 0.0;	B(1,3) = (T) 1.0;
		B(2,0) = q(x,y);  B(2,1) = c(x,y);  B(2,2) = d(x,y);  B(2,3) = em(x,y);
		B(3,0) = C0(x,y); B(3,1) = C1(x,y); B(3,2) = C2(x,y); B(3,3) = C3(x,y);
		// std::cout << "C0: " << C0(x,y) << '\n';
		// std::cout << "C1: " << C1(x,y) << '\n';
		// std::cout << "C2: " << C2(x,y) << '\n';
		// std::cout << "C3: " << C3(x,y) << '\n';
		return true;
	}
} // namespace kv
#endif
