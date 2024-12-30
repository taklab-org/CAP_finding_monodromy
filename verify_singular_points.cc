/*
 * verify_singular_points.cc:
 * This code provides the explicit position of sigular points that are intersection of singular locus and genelic line on C^2.
 *
 * Complile options are as follows (for instance):
 * c++ -I.. -O3 -DNDEBUG verify_singular_points.cc
 * c++ -I.. -O3 -DNDEBUG -DKV_USE_AVX512 -mavx512f -DKV_USE_TPFMA verify_singular_points.cc
 * 
 * written by A. Takayasu	Dec. 27 2024
 */
// 
#include <kv/allsol.hpp>
#include <kv/complex.hpp>
namespace ub = boost::numeric::ublas;
// const int N = 2; // # of complex variables
// Base point: (lam0, mu0)
const kv::interval<double> lam0(pow(2, -10)), mu0(pow(2, -10));
//
struct complex_Func
{ // T is exspected to kv::complex< kv::interval<double> >
	template <class T>
	ub::vector<T> operator()(const ub::vector<T> &x)
	{
		ub::vector<T> y(2);
		T x0, y0;
		x0 = mu0 / pow(lam0 - .25, 2);
		y0 = mu0 / pow(lam0 - .25, 3);
		y(0) = x(0) * x(1) * (4. * x(0) + x(1)) * (pow(36. * x(0) + 27. * x(1) / 2. + 1., 2) - pow(1. - 12. * x(0), 3));
		y(1) = 2. * (x(0) - x0) + (x(1) - y0);
		return y;
	}
};
//
struct Func
{ // T is exspected to kv::interval<double>
	template <class T>
	ub::vector<T> operator()(const ub::vector<T> &x)
	{
		int n = 2; // # of complex variables
		ub::vector<T> y(2 * n);
		complex_Func f;
		ub::vector<kv::complex<T>> x_complex, y_complex;
		x_complex.resize(n);
		// std::cout << x_complex << '\n';
		for (int i = 0; i < n; i++)
		{
			// x_complex(i) = kv::complex<T>(x(i),x(i+n));
			x_complex(i).real() = x(i);
			x_complex(i).imag() = x(i + n);
		}
		// std::cout << "hi!hi!" << '\n';
		y_complex = f(x_complex);
		for (int i = 0; i < n; i++)
		{
			y(i) = y_complex(i).real();
			y(i + n) = y_complex(i).imag();
		}
		return y;
	}
};
int main()
{
	int n = 2; // # of complex variables
	ub::vector<kv::interval<double>> x;
	std::cout.precision(17);
	x.resize(2 * n);
	for (int i = 0; i < 2 * n; i++)
	{
		x(i) = kv::interval<double>("-1000", "1000");
	}
	std::list<ub::vector<kv::interval<double>>> result;
	std::list<ub::vector<kv::interval<double>>>::iterator p;

	result = allsol(Func(), x, 0);

	ub::vector<kv::complex<kv::interval<double>>> singular_points;
	singular_points.resize(n);
	p = result.begin();
	while (p != result.end())
	{
		x = *p;
		for (int i = 0; i < n; i++)
		{
			singular_points(i).real() = x(i);
			singular_points(i).imag() = x(i + n);
		}
		std::cout << "Results: " << singular_points << "\n";
		p++;
	}
}
