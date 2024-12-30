/*
 * initial_values.hpp:
 * Computing initial values via fundamental solutions using series formulation
 *
 * written by A. Takayasu	Dec. 28 2024
 */
//
#ifndef INITIAL_VALUES_HPP
#define INITIAL_VALUES_HPP
//
#include <kv/hwround.hpp>
#include <kv/gamma.hpp>
//
// typedef double dd;
// typedef kv::dd dd;
// typedef kv::interval<dd> intval;
//
// const int N_trunc = 41; // Truncation number
const int N_trunc = 41; // Truncation number
const intval iN(N_trunc);
const kv::complex<intval> im(kv::complex<intval>::i());
const intval ipi(kv::constants<intval>::pi());
const intval ie(kv::constants<intval>::e());
const intval ilog2(kv::constants<intval>::ln2());
const intval abslogmu(10. * ilog2);
//
/////////////////////////// coefficients a, b, c, d ///////////////////////////
// Factorial function for positive integers
intval int_factorial(const int n)
{
	dd p_roundup = (dd)1.0, p_rounddown = (dd)1.0;
	dd p_up, p_down;
	kv::hwround::roundup();
	for (int i = 1; i <= n; i++)
	{
		p_roundup *= (dd)i;
	}
	kv::hwround::rounddown();
	for (int i = 1; i <= n; i++)
	{
		p_rounddown *= (dd)i;
	}
	kv::hwround::roundnear();
	p_up = floor(p_roundup);
	p_down = ceil(p_rounddown);
	if (p_up == p_down)
	{
		return (intval)p_up;
	}
	else
	{
		return intval(p_rounddown, p_roundup);
	}
}
// a_{l,m}
kv::complex<intval> int_alm(const int l, const int m)
{
	kv::complex<intval> a;
	a = int_factorial(2 * l + 4 * m) / int_factorial(l + m) / int_factorial(l) / pow(int_factorial(m), 3);
	return a;
}
// b_{l,m}
kv::complex<intval> int_blm(const int l, const int m)
{
	kv::complex<intval> b;
	intval b4(0.0), b3(0.0);
	for (int j = l + m + 1; j <= 2 * l + 4 * m; j++)
	{
		b4 += 4.0 / (intval)j;
	}
	for (int j = m + 1; j <= l + m; j++)
	{
		b3 += 3.0 / (intval)j;
	}
	b = b4 + b3;
	return b;
}
// c_{l,m}
kv::complex<intval> int_clm(const int l, const int m)
{
	kv::complex<intval> c;
	intval c1(0.0), c2(0.0), c3(0.0);
	for (int j = l + m + 1; j <= 2 * l + 4 * m; j++)
	{
		c1 += 16.0 / pow((intval)j, 2);
	}
	for (int j = m + 1; j <= l + m; j++)
	{
		c2 += 15.0 / pow((intval)j, 2);
	}
	for (int j = 1; j <= m; j++)
	{
		c3 += 12.0 / pow((intval)j, 2);
	}
	c = c1 + c2 + c3;
	return c;
}
// d_{l,m}
kv::complex<intval> int_dlm(const int l, const int m)
{
	kv::complex<intval> d;
	// d = pow(kv::gamma((intval)l+m),3)/kv::gamma((intval)2.*l)/kv::gamma((intval)l+2*m+1)/kv::gamma((intval)m+1.);
	d = pow(int_factorial(l + m - 1), 3) / int_factorial(2 * l - 1) / int_factorial(l + 2 * m) / int_factorial(m);
	if ((l + m) % 2 != 0)
	{
		d = -d;
	}
	return d;
}
// for double input
kv::complex<intval> int_ddlm(const double l, const int m)
{
	kv::complex<intval> d;
	// d = pow((kv::complex<intval>)-1.0,l+m)*pow(kv::gamma((intval)l+m),3)/kv::gamma((intval)2.*l)/kv::gamma((intval)l+2*m+1)/kv::gamma((intval)m+1.);
	d = pow((kv::complex<intval>)-1.0, l + m) * pow(kv::gamma((intval)l + m), 3) / kv::gamma((intval)2. * l) / kv::gamma((intval)l + 2 * m + 1) / int_factorial(m);
	return d;
}
/////////////////////////// truncation errors ///////////////////////////
const intval beta_N(1. + 3. * log(1. + 1. / iN) / (8. * ilog2 + 3. + 3. * log(iN)));
const intval gamm_N(64. * (4. + 1. / iN));
const intval eta1_N((25. / 64.) * (1. + 2. / (2. * iN + 1.)) * sqrt(1. + 1. / iN)); // eta_N^1
const intval eta2_N((25. / 64.) * (1. + 2. / (2. * iN + .5)) * sqrt(1. + 1. / iN)); // eta_N^(1/2)
const intval theta1_N((25. / 64.) * (1. + 1. / (iN + 1.)) * sqrt(1. + 1. / iN));	// theta_N^1
const intval theta2_N((25. / 64.) * (1. + 1. / (iN + .5)) * sqrt(1. + 1. / iN));	// theta_N^(1/2)
const intval iota_N(64. * (4. + 5. / (iN - 1.)));
const intval nu1_N((25. / 64.) * (1. + 2. / (2. * iN + 1.)) * (1. + 2. / (2. * iN)) * sqrt(1. + 1. / iN));		// nu_N^1
const intval nu2_N((25. / 64.) * (1. + 2. / (2. * iN + .5)) * (1. + 2. / (2. * iN - .5)) * sqrt(1. + 1. / iN)); // nu_N^(1/2)
const intval xi1_N((25. / 64.) * (1. + 2. / (2. * iN + 1.)) * (1. + 1. / (iN + 1.)) * sqrt(1. + 1. / iN));		// xi_N^1
const intval xi2_N((25. / 64.) * (1. + 2. / (2. * iN + .5)) * (1. + 1. / (iN + .5)) * sqrt(1. + 1. / iN));		// xi_N^(1/2)
const intval sig1_N((25. / 64.) * (1. + 1. / (iN + 1.)) * (1. + 1. / (iN + 2.)) * sqrt(1. + 1. / iN));			// sigma_N^1
const intval sig2_N((25. / 64.) * (1. + 1. / (iN + .5)) * (1. + 1. / (iN + 1.5)) * sqrt(1. + 1. / iN));			// sigma_N^(1/2)
//
// Suppose (lam0, mu0) = (2^-10, 2^-10)
// Constants:
const intval lam0_mu0(pow(2, -9));
const intval one_plus_lam0(1.0 + (intval)pow(2, -10));
const intval inv_lam0mu0(pow(2, 10));
const intval inv_lam0mu02(pow(2, 20));
const intval N_factorial(int_factorial(N_trunc + 1));
const intval Nm1_factorial(int_factorial(N_trunc));
const intval Nm2_factorial(int_factorial(N_trunc - 1));
const intval FourN_factorial(int_factorial(4 * (N_trunc + 1)));
const intval C1(FourN_factorial / N_factorial / N_factorial / N_factorial / N_factorial);
const intval C2(8. * ilog2 + 3. + 3. * log(iN + 1));
const intval C3(FourN_factorial / N_factorial / N_factorial / N_factorial / Nm1_factorial);
const intval C4(FourN_factorial / N_factorial / N_factorial / N_factorial / Nm2_factorial);
const intval C5(8. * pow(ipi, 2) / 3.);
const intval C6(pow(ie, 3) / (2. * sqrt(2.) * ipi));
// Truncation errors:
const intval delta1(C1 *pow(lam0_mu0, N_trunc + 1) / (1.0 - 256. * lam0_mu0));
const intval delta2(C1 *C2 *pow(lam0_mu0, N_trunc + 1) / (1.0 - 256. * beta_N * lam0_mu0));
const intval delta3(C1 *pow(C2, 2) * pow(lam0_mu0, N_trunc + 1) / (1.0 - 256. * pow(beta_N, 2) * lam0_mu0));
const intval delta4(C5 *C1 *pow(lam0_mu0, N_trunc + 1) / (1.0 - 256. * lam0_mu0));
const intval delta51(C6 *sqrt(iN + 1) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - (25. / 64.) * sqrt(1. + 1. / iN) * one_plus_lam0)); // delta_5^1
const intval delta52(delta51);																													   // dekta_5^(1/2)
const intval delta6(C3 *pow(lam0_mu0, N_trunc) / (1.0 - gamm_N * lam0_mu0));
const intval delta7(C3 *C2 *pow(lam0_mu0, N_trunc) / (1.0 - beta_N * gamm_N * lam0_mu0));
const intval delta8(C3 *pow(C2, 2) * pow(lam0_mu0, N_trunc) / (1.0 - pow(beta_N, 2) * gamm_N * lam0_mu0));
const intval delta9(C5 *C3 *pow(lam0_mu0, N_trunc) / (1.0 - gamm_N * lam0_mu0));
const intval delta101(C6 *inv_lam0mu0 *sqrt(iN + 1) * (2 * iN + 3.) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - eta1_N * one_plus_lam0));	// delta_10^1
const intval delta102(C6 *inv_lam0mu0 *sqrt(iN + 1) * (2 * iN + 2.5) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - eta2_N * one_plus_lam0)); // delta_10^(1/2)
const intval delta111(C6 *inv_lam0mu0 *sqrt(iN + 1) * (iN + 2.) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - theta1_N * one_plus_lam0));	// delta_10^1
const intval delta112(C6 *inv_lam0mu0 *sqrt(iN + 1) * (iN + 1.5) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - theta2_N * one_plus_lam0));	// delta_10^(1/2)
const intval delta12(C4 *pow(lam0_mu0, N_trunc - 1) / (1.0 - iota_N * lam0_mu0));
const intval delta13(C4 *C2 *pow(lam0_mu0, N_trunc - 1) / (1.0 - beta_N * iota_N * lam0_mu0));
const intval delta14(C4 *pow(C2, 2) * pow(lam0_mu0, N_trunc - 1) / (1.0 - pow(beta_N, 2) * iota_N * lam0_mu0));
const intval delta15(C5 *C4 *pow(lam0_mu0, N_trunc - 1) / (1.0 - iota_N * lam0_mu0));
const intval delta161(C6 *inv_lam0mu02 *sqrt(iN + 1) * (2 * iN + 3.) * (2 * iN + 2.) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - nu1_N * one_plus_lam0));
const intval delta162(C6 *inv_lam0mu02 *sqrt(iN + 1) * (2 * iN + 2.5) * (2 * iN + 1.5) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - nu2_N * one_plus_lam0));
const intval delta171(C6 *inv_lam0mu02 *sqrt(iN + 1) * (2 * iN + 3.) * (iN + 2.) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - xi1_N * one_plus_lam0));
const intval delta172(C6 *inv_lam0mu02 *sqrt(iN + 1) * (2 * iN + 2.5) * (iN + 1.5) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - xi2_N * one_plus_lam0));
const intval delta181(C6 *inv_lam0mu02 *sqrt(iN + 1) * (iN + 2.) * (iN + 3.) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - sig1_N * one_plus_lam0));
const intval delta182(C6 *inv_lam0mu02 *sqrt(iN + 1) * (iN + 1.5) * (iN + 2.5) * pow((25. / 64.) * one_plus_lam0, N_trunc + 1) / (1.0 - sig2_N * one_plus_lam0));
//
/////////////////////////// phi1 ///////////////////////////
namespace kv
{
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1(const T lam, const T mu)
	{
		T p;
		p = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p += int_alm(l, m) * pow(lam, l) * pow(mu, m);
			}
		}
		return p + delta1.upper() * intval(-1, 1);
	}
	//----------- First derivatives -----------
	// phi1_lam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_lam(const T lam, const T mu)
	{
		T p_lam;
		p_lam = 0.0;
		for (int l = 1; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lam += l * int_alm(l, m) * pow(lam, l - 1) * pow(mu, m);
			}
		}
		return p_lam + delta6.upper() * intval(-1, 1);
	}
	// phi1_mu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_mu(const T lam, const T mu)
	{
		T p_mu;
		p_mu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 1; m <= N_trunc - l; m++)
			{
				p_mu += m * int_alm(l, m) * pow(lam, l) * pow(mu, m - 1);
			}
		}
		return p_mu + delta6.upper() * intval(-1, 1);
	}
	// lam_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_lam_x(const T lam, const T mu)
	{
		T y;
		y = mu / pow(lam - 0.25, 3);
		return 1.0 / y;
	}
	// mu_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_mu_x(const T lam, const T mu)
	{
		T x, y;
		x = mu / pow(lam - 0.25, 2);
		y = mu / pow(lam - 0.25, 3);
		return 3.0 * pow(x, 2) / pow(y, 2);
	}
	// lam_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_lam_y(const T lam, const T mu)
	{
		T x, y;
		x = mu / pow(lam - 0.25, 2);
		y = mu / pow(lam - 0.25, 3);
		return -x / pow(y, 2);
	}
	// mu_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_mu_y(const T lam, const T mu)
	{
		T x, y;
		x = mu / pow(lam - 0.25, 2);
		y = mu / pow(lam - 0.25, 3);
		return -2.0 * pow(x, 3) / pow(y, 3);
	}
	// phi1_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_x(const T lam, const T mu)
	{
		return int_phi1_lam(lam, mu) * int_lam_x(lam, mu) + int_phi1_mu(lam, mu) * int_mu_x(lam, mu);
	}
	// phi1_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_y(const T lam, const T mu)
	{
		return int_phi1_lam(lam, mu) * int_lam_y(lam, mu) + int_phi1_mu(lam, mu) * int_mu_y(lam, mu);
	}
	//----------- Second derivative -----------
	// lam_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_lam_xy(const T lam, const T mu)
	{
		T y;
		y = mu / pow(lam - 0.25, 3);
		return -1.0 / pow(y, 2);
	}
	// mu_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_mu_xy(const T lam, const T mu)
	{
		T x, y;
		x = mu / pow(lam - 0.25, 2);
		y = mu / pow(lam - 0.25, 3);
		return -6.0 * pow(x, 2) / pow(y, 3);
	}
	// phi1_lamlam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_lamlam(const T lam, const T mu)
	{
		T p_lamlam;
		p_lamlam = 0.0;
		for (int l = 2; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lamlam += l * (l - 1) * int_alm(l, m) * pow(lam, l - 2) * pow(mu, m);
			}
		}
		return p_lamlam + delta12.upper() * intval(-1, 1);
	}
	// phi1_lammu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_lammu(const T lam, const T mu)
	{
		T p_lammu;
		p_lammu = 0.0;
		for (int l = 1; l <= N_trunc; l++)
		{
			for (int m = 1; m <= N_trunc - l; m++)
			{
				p_lammu += l * m * int_alm(l, m) * pow(lam, l - 1) * pow(mu, m - 1);
			}
		}
		return p_lammu + delta12.upper() * intval(-1, 1);
	}
	// phi1_mumu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_mumu(const T lam, const T mu)
	{
		T p_mumu;
		p_mumu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 2; m <= N_trunc - l; m++)
			{
				p_mumu += m * (m - 1) * int_alm(l, m) * pow(lam, l) * pow(mu, m - 2);
			}
		}
		return p_mumu + delta12.upper() * intval(-1, 1);
	}
	// phi1_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi1_xy(const T lam, const T mu)
	{
		T phi_lam_x = int_phi1_lamlam(lam, mu) * int_lam_x(lam, mu) + int_phi1_lammu(lam, mu) * int_mu_x(lam, mu);
		T phi_mu_x = int_phi1_lammu(lam, mu) * int_lam_x(lam, mu) + int_phi1_mumu(lam, mu) * int_mu_x(lam, mu);
		return (phi_lam_x)*int_lam_y(lam, mu) + int_phi1_lam(lam, mu) * int_lam_xy(lam, mu) + (phi_mu_x)*int_mu_y(lam, mu) + int_phi1_mu(lam, mu) * int_mu_xy(lam, mu);
	}
	//
	/////////////////////////// phi2 ///////////////////////////
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2(const T lam, const T mu)
	{
		T p;
		p = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p += int_alm(l, m) * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m)) * pow(lam, l) * pow(mu, m) + 0.5 * int_dlm(l + 1, m) * pow(lam, l + 2 * m + 1) / pow(mu, l + m + 1);
			}
		}
		intval err((pow(abslogmu, 2) * delta1 + 2 * abslogmu * delta2 + delta3 + delta4 + 0.5 * delta51) / (2 * pow(ipi, 2)));
		return p / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	//----------- First derivatives -----------
	// phi2_lam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_lam(const T lam, const T mu)
	{
		T p_lam;
		p_lam = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lam += l * int_alm(l, m) * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m)) * pow(lam, l - 1) * pow(mu, m) + 0.5 * int_dlm(l + 1, m) * (l + 2 * m + 1) * pow(lam, l + 2 * m) / pow(mu, l + m + 1);
			}
		}
		intval err((pow(abslogmu, 2) * delta6 + 2 * abslogmu * delta7 + delta8 + delta9 + 0.5 * delta101) / (2 * pow(ipi, 2)));
		return p_lam / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	// phi2_mu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_mu(const T lam, const T mu)
	{
		T p_mu;
		p_mu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mu += int_alm(l, m) * (2 * (log(mu) + int_blm(l, m)) + m * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m))) * pow(lam, l) * pow(mu, m - 1) - 0.5 * int_dlm(l + 1, m) * (l + m + 1) * pow(lam, l + 2 * m + 1) / pow(mu, l + m + 2);
			}
		}
		intval err((2 * abslogmu / abs(mu) * delta1 + 2.0 / abs(mu) * delta2 + pow(abslogmu, 2) * delta6 + 2 * abslogmu * delta7 + delta8 + delta9 + 0.5 * delta111) / (2 * pow(ipi, 2)));
		return p_mu / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	// phi2_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_x(const T lam, const T mu)
	{
		return int_phi2_lam(lam, mu) * int_lam_x(lam, mu) + int_phi2_mu(lam, mu) * int_mu_x(lam, mu);
	}
	// phi2_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_y(const T lam, const T mu)
	{
		return int_phi2_lam(lam, mu) * int_lam_y(lam, mu) + int_phi2_mu(lam, mu) * int_mu_y(lam, mu);
	}
	//----------- Second derivative -----------
	// phi2_lamlam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_lamlam(const T lam, const T mu)
	{
		T p_lamlam;
		p_lamlam = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lamlam += l * (l - 1) * int_alm(l, m) * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m)) * pow(lam, l - 2) * pow(mu, m) + 0.5 * int_dlm(l + 1, m) * (l + 2 * m + 1) * (l + 2 * m) * pow(lam, l + 2 * m - 1) / pow(mu, l + m + 1);
			}
		}
		intval err((pow(abslogmu, 2) * delta12 + 2 * abslogmu * delta13 + delta14 + delta15 + 0.5 * delta161) / (2 * pow(ipi, 2)));
		return p_lamlam / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	// phi2_lammu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_lammu(const T lam, const T mu)
	{
		T p_lammu;
		p_lammu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lammu += l * int_alm(l, m) * (2 * (log(mu) + int_blm(l, m)) + m * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m))) * pow(lam, l - 1) * pow(mu, m - 1) - 0.5 * int_dlm(l + 1, m) * (l + m + 1) * (l + 2 * m + 1) * pow(lam, l + 2 * m) / pow(mu, l + m + 2);
			}
		}
		intval err((2 * abslogmu / abs(mu) * delta6 + 2.0 / abs(mu) * delta7 + pow(abslogmu, 2) * delta12 + 2 * abslogmu * delta13 + delta14 + delta15 + 0.5 * delta171) / (2 * pow(ipi, 2)));
		return p_lammu / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	// phi2_mumu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_mumu(const T lam, const T mu)
	{
		T p_mumu;
		p_mumu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mumu += int_alm(l, m) * (2 + (4 * m - 2) * (log(mu) + int_blm(l, m)) + m * (m - 1) * (pow(log(mu) + int_blm(l, m), 2) - int_clm(l, m))) * pow(lam, l) * pow(mu, m - 2) + 0.5 * int_dlm(l + 1, m) * (l + m + 1) * (l + m + 2) * pow(lam, l + 2 * m + 1) / pow(mu, l + m + 3);
			}
		}
		intval err(2 / pow(abs(mu), 2) * delta1 + (4 * abslogmu / abs(mu) * delta6 + 4.0 / abs(mu) * delta7 + pow(abslogmu, 2) * delta12 + 2 * abslogmu * delta13 + delta14 + delta15 + 0.5 * delta181) / (2 * pow(ipi, 2)));
		return p_mumu / (2 * pow(ipi, 2)) + err.upper() * intval(-1, 1);
	}
	// phi2_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi2_xy(const T lam, const T mu)
	{
		T phi_lam_x = int_phi2_lamlam(lam, mu) * int_lam_x(lam, mu) + int_phi2_lammu(lam, mu) * int_mu_x(lam, mu);
		T phi_mu_x = int_phi2_lammu(lam, mu) * int_lam_x(lam, mu) + int_phi2_mumu(lam, mu) * int_mu_x(lam, mu);
		return (phi_lam_x)*int_lam_y(lam, mu) + int_phi2_lam(lam, mu) * int_lam_xy(lam, mu) + (phi_mu_x)*int_mu_y(lam, mu) + int_phi2_mu(lam, mu) * int_mu_xy(lam, mu);
	}
	//
	/////////////////////////// phi3 ///////////////////////////
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3(const T lam, const T mu)
	{
		T p;
		p = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p += int_ddlm(l + 0.5, m) * pow(lam, l + 2 * m + 0.5) / pow(mu, l + m + 0.5);
			}
		}
		intval err(delta52 / (4 * pow(ipi, 2)));
		return p / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	//----------- First derivatives -----------
	// phi3_lam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_lam(const T lam, const T mu)
	{
		T p_lam;
		p_lam = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lam += int_ddlm(l + 0.5, m) * (l + 2 * m + 0.5) * pow(lam, l + 2 * m - 0.5) / pow(mu, l + m + 0.5);
			}
		}
		intval err(delta102 / (4 * pow(ipi, 2)));
		return p_lam / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	// phi3_mu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_mu(const T lam, const T mu)
	{
		T p_mu;
		p_mu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mu += int_ddlm(l + 0.5, m) * (l + m + 0.5) * pow(lam, l + 2 * m + 0.5) / pow(mu, l + m + 1.5);
			}
		}
		intval err(delta112 / (4 * pow(ipi, 2)));
		return -p_mu / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	// phi3_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_x(const T lam, const T mu)
	{
		return int_phi3_lam(lam, mu) * int_lam_x(lam, mu) + int_phi3_mu(lam, mu) * int_mu_x(lam, mu);
	}
	// phi3_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_y(const T lam, const T mu)
	{
		return int_phi3_lam(lam, mu) * int_lam_y(lam, mu) + int_phi3_mu(lam, mu) * int_mu_y(lam, mu);
	}
	//----------- Second derivative -----------
	// phi3_lamlam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_lamlam(const T lam, const T mu)
	{
		T p_lamlam;
		p_lamlam = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lamlam += int_ddlm(l + 0.5, m) * (l + 2 * m + 0.5) * (l + 2 * m - 0.5) * pow(lam, l + 2 * m - 1.5) / pow(mu, l + m + 0.5);
			}
		}
		intval err(delta162 / (4 * pow(ipi, 2)));
		return p_lamlam / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	// phi3_lammu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_lammu(const T lam, const T mu)
	{
		T p_lammu;
		p_lammu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lammu += int_ddlm(l + 0.5, m) * (l + 2 * m + 0.5) * (l + m + 0.5) * pow(lam, l + 2 * m - 0.5) / pow(mu, l + m + 1.5);
			}
		}
		intval err(delta172 / (4 * pow(ipi, 2)));
		return -p_lammu / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	// phi3_mumu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_mumu(const T lam, const T mu)
	{
		T p_mumu;
		p_mumu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mumu += int_ddlm(l + 0.5, m) * (l + m + 0.5) * (l + m + 1.5) * pow(lam, l + 2 * m + 0.5) / pow(mu, l + m + 2.5);
			}
		}
		intval err(delta182 / (4 * pow(ipi, 2)));
		return p_mumu / (4 * pow(ipi, 2)) + err.upper() * im;
	}
	// phi3_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi3_xy(const T lam, const T mu)
	{
		T lamx = int_lam_x(lam, mu);
		T mux = int_mu_x(lam, mu);
		T phi3_lammu = int_phi3_lammu(lam, mu);
		T phi_lam_x = int_phi3_lamlam(lam, mu) * lamx + phi3_lammu * mux;
		T phi_mu_x = phi3_lammu * lamx + int_phi3_mumu(lam, mu) * mux;
		return (phi_lam_x)*int_lam_y(lam, mu) + int_phi3_lam(lam, mu) * int_lam_xy(lam, mu) + (phi_mu_x)*int_mu_y(lam, mu) + int_phi3_mu(lam, mu) * int_mu_xy(lam, mu);
	}
	//
	/////////////////////////// phi4 ///////////////////////////
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4(const T lam, const T mu)
	{
		T p;
		p = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p += int_alm(l, m) * (log(mu) + int_blm(l, m)) * pow(lam, l) * pow(mu, m);
			}
		}
		intval err((abslogmu * delta1 + delta2) / (2 * ipi));
		return p / (2 * ipi * im) + err.upper() * im;
	}
	//----------- First derivatives -----------
	// phi4_lam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_lam(const T lam, const T mu)
	{
		T p_lam;
		p_lam = 0.0;
		for (int l = 1; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lam += l * int_alm(l, m) * (log(mu) + int_blm(l, m)) * pow(lam, l - 1) * pow(mu, m);
			}
		}
		intval err((abslogmu * delta6 + delta7) / (2 * ipi));
		return p_lam / (2 * ipi * im) + err.upper() * im;
	}
	// phi4_mu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_mu(const T lam, const T mu)
	{
		T p_mu;
		p_mu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mu += int_alm(l, m) * (1.0 + m * (log(mu) + int_blm(l, m))) * pow(lam, l) * pow(mu, m - 1);
			}
		}
		intval err((1.0 / abs(mu) * delta1 + abslogmu * delta6 + delta7) / (2 * ipi));
		return p_mu / (2 * ipi * im) + err.upper() * im;
	}
	// phi4_x
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_x(const T lam, const T mu)
	{
		return int_phi4_lam(lam, mu) * int_lam_x(lam, mu) + int_phi4_mu(lam, mu) * int_mu_x(lam, mu);
	}
	// phi4_y
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_y(const T lam, const T mu)
	{
		return int_phi4_lam(lam, mu) * int_lam_y(lam, mu) + int_phi4_mu(lam, mu) * int_mu_y(lam, mu);
	}
	//----------- Second derivative -----------
	// phi4_lamlam
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_lamlam(const T lam, const T mu)
	{
		T p_lamlam;
		p_lamlam = 0.0;
		for (int l = 2; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lamlam += l * (l - 1) * int_alm(l, m) * (log(mu) + int_blm(l, m)) * pow(lam, l - 2) * pow(mu, m);
			}
		}
		intval err((abslogmu * delta12 + delta13) / (2 * ipi));
		return p_lamlam / (2 * ipi * im) + err.upper() * im;
	}
	// phi4_lammu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_lammu(const T lam, const T mu)
	{
		T p_lammu;
		p_lammu = 0.0;
		for (int l = 1; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_lammu += l * int_alm(l, m) * (1.0 + m * (log(mu) + int_blm(l, m))) * pow(lam, l - 1) * pow(mu, m - 1);
			}
		}
		intval err((1.0 / abs(mu) * delta6 + abslogmu * delta12 + delta13) / (2 * ipi));
		return p_lammu / (2 * ipi * im) + err.upper() * im;
	}
	// phi4_mumu
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_mumu(const T lam, const T mu)
	{
		T p_mumu;
		p_mumu = 0.0;
		for (int l = 0; l <= N_trunc; l++)
		{
			for (int m = 0; m <= N_trunc - l; m++)
			{
				p_mumu += int_alm(l, m) * ((m - 1) + m * (1 + (m - 1) * (log(mu) + int_blm(l, m)))) * pow(lam, l) * pow(mu, m - 2);
			}
		}
		intval err((2.0 / abs(mu) * delta6 + abslogmu * delta12 + delta13) / (2 * ipi));
		return p_mumu / (2 * ipi * im) + err.upper() * im;
	}
	// phi4_xy
	// T is expected to kv::complex< kv::interval<double> >
	template <class T>
	T int_phi4_xy(const T lam, const T mu)
	{
		T phi_lam_x = int_phi4_lamlam(lam, mu) * int_lam_x(lam, mu) + int_phi4_lammu(lam, mu) * int_mu_x(lam, mu);
		T phi_mu_x = int_phi4_lammu(lam, mu) * int_lam_x(lam, mu) + int_phi4_mumu(lam, mu) * int_mu_x(lam, mu);
		return (phi_lam_x)*int_lam_y(lam, mu) + int_phi4_lam(lam, mu) * int_lam_xy(lam, mu) + (phi_mu_x)*int_mu_y(lam, mu) + int_phi4_mu(lam, mu) * int_mu_xy(lam, mu);
	}
} // namespace kv
#endif
