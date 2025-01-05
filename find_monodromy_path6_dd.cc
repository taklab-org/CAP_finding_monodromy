/*
 * find_monodromy_path6_dd.cc:
 * This code provides the result of rigorous analytic continuation from the base point using interval arithmeri with DD precision.
 *
 * Complile options are as follows (for instance):
 * c++ -I.. -O3 -DNDEBUG find_monodromy_path6_dd.cc -fopenmp
 * c++ -I.. -O3 -DNDEBUG -DKV_USE_AVX512 -mavx512f -DKV_USE_TPFMA find_monodromy_path6_dd.cc -fopenmp
 * 
 * written by A. Takayasu	Dec. 29 2024
 */
// 
#include <kv/ode-maffine.hpp>
#include <kv/complex.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/vleq.hpp>
#include <omp.h>
#include "variable_coefficients.hpp"
// 
namespace ub = boost::numeric::ublas;
// typedef double dd;
typedef kv::dd dd;
typedef kv::interval<dd> intval;
// 
#include "initial_values.hpp"
// 
const int N = 4; // dimention of ODE (# of complex variables)
const int order_of_taylor = 15; // maximal order of power series
const kv::complex<intval> lam0(pow(2,-10)), mu0(pow(2,-10));
// 
// Complex vector fields of the Maurer-Cartan form on a (semi-)circle
template <class TT> struct odefunc_anticlockwise {// T is expected to kv::complex<intval>
  TT cx, cy, r;
  odefunc_anticlockwise(TT cx, TT cy, TT r): cx(cx), cy(cy), r(r) {}

  template <class T> ub::vector<T> operator() (const ub::vector<T>& phi, const T& t) {
    int n = N;
    T im = (T) kv::complex<intval>::i();
    T ipi = (T) kv::constants<intval>::pi();
    T expipi = exp(im*ipi*t);
    T x, y, dx, dy;
    // path: centered at c = (cx,cy) with radious r
    x = (T) cx + (T)r * expipi;
    y = (T) cy - 2. * (T) r * expipi;
    dx = im * ipi * (T) r * expipi;
    dy = -2. * dx;
    // general
    ub::vector<T> dphi;
    dphi.resize(n);
    ub::matrix<T> A, B, C;
    A.resize(n,n);
    B.resize(n,n);
    matA(x,y,A);
    matB(x,y,B);
    C = dx*A + dy*B;
    dphi = prod(C,phi);
    return dphi;
  }
};
// Wrapper of complex function "odefunc_anticlockwise" for odelong_maffine
template <class TT> struct vec_odefunc_anticlockwise {//  T is expected to kv::interval<double>
  TT cx, cy, r;
  vec_odefunc_anticlockwise(TT cx, TT cy, TT r): cx(cx), cy(cy), r(r) {}

  template <class T> ub::vector<T> operator() (const ub::vector<T>& x, const T& t) {
    // odefunc_anticlockwise<TT> f_complex(cx,cy,r,wgt);
    odefunc_anticlockwise<TT> f_complex(cx,cy,r);
    int n = N; // # of complex variables
    ub::vector< kv::complex<T> > x_complex, y_complex;
    x_complex.resize(n);
    // std::cout << x_complex << '\n';
    for (int i = 0; i < n; i++) {
      // x_complex(i) = kv::complex<T>(x(i),x(i+n));
      x_complex(i).real() = x(i);
      x_complex(i).imag() = x(i+n);
    }
    y_complex = f_complex(x_complex, (kv::complex<T>)t);
    ub::vector<T> y(2*n);
    for (int i = 0; i < n; i++) {
      y(i) = y_complex(i).real();
      y(i+n) = y_complex(i).imag();
    }
    return y;
  }
};
// Complex vector fields of the Maurer-Cartan form on a segment
template <class TT> struct odefunc_segment {// T is expected to kv::complex<intval>
  TT ax, ay, bx, by;
  odefunc_segment(TT ax, TT ay, TT bx, TT by): ax(ax), ay(ay), bx(bx), by(by) {}

  template <class T> ub::vector<T> operator() (const ub::vector<T>& phi, const T& t) {
    int n = N;
    // T im = (T) kv::complex<intval>::i();
    // T ipi = (T) kv::constants<intval>::pi();
    // T expipi = exp(-im*ipi*t);
    T x, y, dx, dy;
    // path from a to b
    x = (1-t) * (T) ax + t * (T) bx;
    y = (1-t) * (T) ay + t * (T) by;
    dx = bx - ax;
    dy = by - ay;
    // x = (T) cx + (T) r * expipi;
    // y = (T) cy - 2. * (T) r * expipi;
    // dx = - im * ipi * (T) r * expipi;
    // dy = -2. * dx;
    // general
    ub::vector<T> dphi;
    dphi.resize(n);
    ub::matrix<T> A, B, C;
    A.resize(n,n);
    B.resize(n,n);
    matA(x,y,A);
    matB(x,y,B);
    C = dx * A + dy * B;
    dphi = prod(C,phi);
    return dphi;
  }
};
// Wrapper of complex function "odefunc_segment" for odelong_maffine
template <class TT> struct vec_odefunc_segment {//  T is expected to kv::interval<double>
  TT ax, ay, bx, by;
  vec_odefunc_segment(TT ax, TT ay, TT bx, TT by): ax(ax), ay(ay), bx(bx), by(by) {}

  template <class T> ub::vector<T> operator() (const ub::vector<T>& x, const T& t) {
    odefunc_segment<TT> f_complex(ax,ay,bx,by);
    // 
    int n = N; // # of complex variables
    ub::vector< kv::complex<T> > x_complex, y_complex;
    x_complex.resize(n);
    for (int i = 0; i < n; i++) {
      x_complex(i).real() = x(i);
      x_complex(i).imag() = x(i+n);
    }
    y_complex = f_complex(x_complex, (kv::complex<T>)t);
    ub::vector<T> y(2*n);
    for (int i = 0; i < n; i++) {
      y(i) = y_complex(i).real();
      y(i+n) = y_complex(i).imag();
    }
    return y;
  }
};
//
int main()
{
std::cout << "Max threads: " << omp_get_max_threads() << '\n';
std::cout.precision(17);
int n = N; // # of complex variables
int r;
/////// Initial value via fundamental solution phi1, phi2, phi3, phi4
ub::vector< kv::complex<intval> > phi1, phi1_init;
ub::vector< kv::complex<intval> > phi2, phi2_init;
ub::vector< kv::complex<intval> > phi3, phi3_init;
ub::vector< kv::complex<intval> > phi4, phi4_init;
phi1.resize(n);
phi2.resize(n);
phi3.resize(n);
phi4.resize(n);
//
/////// Define the path gamma_{6,j} (j=1,2,3)
kv::complex<intval> x0, y0, x4, y4, x6, y6;
kv::complex<intval> c6x, c6y, r6;
x0 = mu0/(pow(lam0-0.25,2)); y0 = mu0/(pow(lam0-0.25,3)); // Base point
// x1 = -0.015871301384836914; y1 = 0; // p1 singular point
// x2 = 0; y2 = -0.031742602769673827; // p2 singular point
x4 = 0.016430419036184044; y4 = -0.064603440842041781; // p4 singular point
x6 = kv::complex<intval>(0.093347290481908138, -0.12249493723489422);
y6 = kv::complex<intval>(-0.21843718373348947, 0.24498987446978776); // p6 singular point
r6 = abs(x6-x4)/2;
c6x = x6-r6;
c6y = y6+2*r6;
// parallel implementation for analytic continuation
intval start1, end1;
intval start2, end2;
intval start3, end3;
intval start4, end4;
ub::vector<intval> phi1_vec;
ub::vector<intval> phi2_vec;
ub::vector<intval> phi3_vec;
ub::vector<intval> phi4_vec;
ub::vector< kv::complex<intval> > phi1_path61, phi1_path62, phi1_path63;
ub::vector< kv::complex<intval> > phi2_path61, phi2_path62, phi2_path63;
ub::vector< kv::complex<intval> > phi3_path61, phi3_path62, phi3_path63;
ub::vector< kv::complex<intval> > phi4_path61, phi4_path62, phi4_path63;
//////////////////////////////////////////////////////////////// parallel start ///////////////////////////////////////////////////
#pragma omp parallel sections
{
///////////////////////////////////////////////////////////////// Computing phi1 ///////////////////////////////////////////////////
#pragma omp section
{
  phi1(0) = int_phi1(lam0,mu0);
  phi1(1) = int_phi1_x(lam0,mu0);
  phi1(2) = int_phi1_y(lam0,mu0);
  phi1(3) = int_phi1_xy(lam0,mu0);
  // std::cout << "phi1: " << phi1 << '\n';
  phi1_init = phi1;
  std::cout << "phi1_init: " << phi1_init << '\n';
  // Initial values for ODEs
  phi1(0) = kv::complex<intval>(1,0);
  phi1(1) = kv::complex<intval>(0,0);
  phi1(2) = kv::complex<intval>(0,0);
  phi1(3) = kv::complex<intval>(0,0);
  //
  /////// Compute the analytic continuation on the path gamma_{6,1}
  phi1_vec.resize(2*n);
  //
  for (int i = 0; i < n; i++) {
    phi1_vec(i) = phi1(i).real();
    phi1_vec(i+n) = phi1(i).imag();
  }
  // std::cout << "phi_vec: " << phi1_vec << '\n';
  vec_odefunc_segment< kv::complex<intval> > f61(x0,y0,c6x,c6y);
  start1 = 0;
  end1 = 1;
  std::cout << "Computing phi1 along the path gamma_{6,1}..." << '\n';
  r = kv::odelong_maffine(
    f61,
    phi1_vec,
    start1,
    end1,
    // kv::ode_param<double>().set_verbose(1).set_order(24)
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi1_path61: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi1_path61: Solution incomplitely calculated until t = " << end1 << ".\n";
  } else if (r == 3) {
    std::cout << "phi1_path61: Solution calculated until t = " << end1 << ".\n";
  } else {
    std::cout << "phi1_path61: Solution completely calculated until t = " << end1 << ".\n";
  }
  //
  phi1_path61.resize(n);
  for (int i = 0; i < n; i++) {
    phi1_path61(i).real() = phi1_vec(i);
    phi1_path61(i).imag() = phi1_vec(i+n);
  }
  std::cout << "phi1_path61: " << phi1_path61 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,2}
  vec_odefunc_anticlockwise< kv::complex<intval> > f62(x6,y6,r6);
  start1 = -1;
  end1 = 1;
  std::cout << "Computing phi1 along the path gamma_{6,2}..." << '\n';
  r = kv::odelong_maffine(
    f62,
    phi1_vec,
    start1,
    end1,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi1_path62: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi1_path62: Solution incomplitely calculated until t = " << end1 << ".\n";
  } else if (r == 3) {
    std::cout << "phi1_path62: Solution calculated until t = " << end1 << ".\n";
  } else {
    std::cout << "phi1_path62: Solution completely calculated until t = " << end1 << ".\n";
  }
  //
  phi1_path62.resize(n);
  for (int i = 0; i < n; i++) {
    phi1_path62(i).real() = phi1_vec(i);
    phi1_path62(i).imag() = phi1_vec(i+n);
  }
  std::cout << "phi1_path62: " << phi1_path62 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,3}
  vec_odefunc_segment< kv::complex<intval> > f63(c6x,c6y,x0,y0);
  start1 = 0;
  end1 = 1;
  std::cout << "Computing phi1 along the path gamma_{6,3}..." << '\n';
  r = kv::odelong_maffine(
    f63,
    phi1_vec,
    start1,
    end1,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi1_path63: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi1_path63: Solution incomplitely calculated until t = " << end1 << ".\n";
  } else if (r == 3) {
    std::cout << "phi1_path63: Solution calculated until t = " << end1 << ".\n";
  } else {
    std::cout << "phi1_path63: Solution completely calculated until t = " << end1 << ".\n";
  }
  //
  phi1_path63.resize(n);
  for (int i = 0; i < n; i++) {
    phi1_path63(i).real() = phi1_vec(i);
    phi1_path63(i).imag() = phi1_vec(i+n);
  }
  std::cout << "phi1_final: " << phi1_path63 << '\n';
  std::cout << "thread = " << omp_get_thread_num() << '\n';
}
///////////////////////////////////////////////////////////////// Computing phi2 ///////////////////////////////////////////////////
#pragma omp section
{
  phi2(0) = int_phi2(lam0,mu0);
  phi2(1) = int_phi2_x(lam0,mu0);
  phi2(2) = int_phi2_y(lam0,mu0);
  phi2(3) = int_phi2_xy(lam0,mu0);
  // std::cout << "phi2: " << phi2 << '\n';
  phi2_init = phi2;
  std::cout << "phi2_init: " << phi2_init << '\n';
  // Initial values for ODEs
  phi2(0) = kv::complex<intval>(0,0);
  phi2(1) = kv::complex<intval>(1,0);
  phi2(2) = kv::complex<intval>(0,0);
  phi2(3) = kv::complex<intval>(0,0);
  //
  /////// Compute the analytic continuation on the path gamma_{6,1}
  phi2_vec.resize(2*n);
  //
  for (int i = 0; i < n; i++) {
    phi2_vec(i) = phi2(i).real();
    phi2_vec(i+n) = phi2(i).imag();
  }
  // std::cout << "phi_vec: " << phi2_vec << '\n';
  vec_odefunc_segment< kv::complex<intval> > f61(x0,y0,c6x,c6y);
  start2 = 0;
  end2 = 1;
  std::cout << "Computing phi2 along the path gamma_{6,1}..." << '\n';
  r = kv::odelong_maffine(
    f61,
    phi2_vec,
    start2,
    end2,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi2_path61: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi2_path61: Solution incomplitely calculated until t = " << end2 << ".\n";
  } else if (r == 3) {
    std::cout << "phi2_path61: Solution calculated until t = " << end2 << ".\n";
  } else {
    std::cout << "phi2_path61: Solution completely calculated until t = " << end2 << ".\n";
  }
  //
  phi2_path61.resize(n);
  for (int i = 0; i < n; i++) {
    phi2_path61(i).real() = phi2_vec(i);
    phi2_path61(i).imag() = phi2_vec(i+n);
  }
  std::cout << "phi2_path61: " << phi2_path61 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,2}
  vec_odefunc_anticlockwise< kv::complex<intval> > f62(x6,y6,r6);
  start2 = -1;
  end2 = 1;
  std::cout << "Computing phi2 along the path gamma_{6,2}..." << '\n';
  r = kv::odelong_maffine(
    f62,
    phi2_vec,
    start2,
    end2,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi2_path62: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi2_path62: Solution incomplitely calculated until t = " << end2 << ".\n";
  } else if (r == 3) {
    std::cout << "phi2_path62: Solution calculated until t = " << end2 << ".\n";
  } else {
    std::cout << "phi2_path62: Solution completely calculated until t = " << end2 << ".\n";
  }
  //
  phi2_path62.resize(n);
  for (int i = 0; i < n; i++) {
    phi2_path62(i).real() = phi2_vec(i);
    phi2_path62(i).imag() = phi2_vec(i+n);
  }
  std::cout << "phi2_path62: " << phi2_path62 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,3}
  vec_odefunc_segment< kv::complex<intval> > f63(c6x,c6y,x0,y0);
  start2 = 0;
  end2 = 1;
  std::cout << "Computing phi2 along the path gamma_{6,3}..." << '\n';
  r = kv::odelong_maffine(
    f63,
    phi2_vec,
    start2,
    end2,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi2_path63: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi2_path63: Solution incomplitely calculated until t = " << end2 << ".\n";
  } else if (r == 3) {
    std::cout << "phi2_path63: Solution calculated until t = " << end2 << ".\n";
  } else {
    std::cout << "phi2_path63: Solution completely calculated until t = " << end2 << ".\n";
  }
  //
  phi2_path63.resize(n);
  for (int i = 0; i < n; i++) {
    phi2_path63(i).real() = phi2_vec(i);
    phi2_path63(i).imag() = phi2_vec(i+n);
  }
  std::cout << "phi2_final: " << phi2_path63 << '\n';
  std::cout << "thread = " << omp_get_thread_num() << '\n';
}
///////////////////////////////////////////////////////////////// Computing phi3 ///////////////////////////////////////////////////
#pragma omp section
{
  phi3(0) = int_phi3(lam0,mu0); // this takes more than 5mins
  phi3(1) = int_phi3_x(lam0,mu0); // this takes more than 5mins
  phi3(2) = int_phi3_y(lam0,mu0); // this takes more than 5mins
  phi3(3) = int_phi3_xy(lam0,mu0); // this takes more than 5mins
  // std::cout << "phi3: " << phi3 << '\n';
  phi3_init = phi3;
  std::cout << "phi3_init: " << phi3_init << '\n';
  // Initial values for ODEs
  phi3(0) = kv::complex<intval>(0,0);
  phi3(1) = kv::complex<intval>(0,0);
  phi3(2) = kv::complex<intval>(1,0);
  phi3(3) = kv::complex<intval>(0,0);
  //
  /////// Compute the analytic continuation on the path gamma_{6,1}
  phi3_vec.resize(2*n);
  //
  for (int i = 0; i < n; i++) {
    phi3_vec(i) = phi3(i).real();
    phi3_vec(i+n) = phi3(i).imag();
  }
  // std::cout << "phi_vec: " << phi3_vec << '\n';
  vec_odefunc_segment< kv::complex<intval> > f61(x0,y0,c6x,c6y);
  start3 = 0;
  end3 = 1;
  std::cout << "Computing phi3 along the path gamma_{6,1}..." << '\n';
  r = kv::odelong_maffine(
    f61,
    phi3_vec,
    start3,
    end3,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi3_path61: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi3_path61: Solution incomplitely calculated until t = " << end3 << ".\n";
  } else if (r == 3) {
    std::cout << "phi3_path61: Solution calculated until t = " << end3 << ".\n";
  } else {
    std::cout << "phi3_path61: Solution completely calculated until t = " << end3 << ".\n";
  }
  //
  phi3_path61.resize(n);
  for (int i = 0; i < n; i++) {
    phi3_path61(i).real() = phi3_vec(i);
    phi3_path61(i).imag() = phi3_vec(i+n);
  }
  std::cout << "phi3_path61: " << phi3_path61 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,2}
  vec_odefunc_anticlockwise< kv::complex<intval> > f62(x6,y6,r6);
  start3 = -1;
  end3 = 1;
  std::cout << "Computing phi3 along the path gamma_{6,2}..." << '\n';
  r = kv::odelong_maffine(
    f62,
    phi3_vec,
    start3,
    end3,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi3_path62: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi3_path62: Solution incomplitely calculated until t = " << end3 << ".\n";
  } else if (r == 3) {
    std::cout << "phi3_path62: Solution calculated until t = " << end3 << ".\n";
  } else {
    std::cout << "phi3_path62: Solution completely calculated until t = " << end3 << ".\n";
  }
  //
  phi3_path62.resize(n);
  for (int i = 0; i < n; i++) {
    phi3_path62(i).real() = phi3_vec(i);
    phi3_path62(i).imag() = phi3_vec(i+n);
  }
  std::cout << "phi3_path62: " << phi3_path62 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,3}
  vec_odefunc_segment< kv::complex<intval> > f63(c6x,c6y,x0,y0);
  start3 = 0;
  end3 = 1;
  std::cout << "Computing phi3 along the path gamma_{6,3}..." << '\n';
  r = kv::odelong_maffine(
    f63,
    phi3_vec,
    start3,
    end3,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi3_path63: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi3_path63: Solution incomplitely calculated until t = " << end3 << ".\n";
  } else if (r == 3) {
    std::cout << "phi3_path63: Solution calculated until t = " << end3 << ".\n";
  } else {
    std::cout << "phi3_path63: Solution completely calculated until t = " << end3 << ".\n";
  }
  //
  phi3_path63.resize(n);
  for (int i = 0; i < n; i++) {
    phi3_path63(i).real() = phi3_vec(i);
    phi3_path63(i).imag() = phi3_vec(i+n);
  }
  std::cout << "phi3_final: " << phi3_path63 << '\n';
  std::cout << "thread = " << omp_get_thread_num() << '\n';
}
///////////////////////////////////////////////////////////////// Computing phi4 ///////////////////////////////////////////////////
#pragma omp section
{
  phi4(0) = int_phi4(lam0,mu0);
  phi4(1) = int_phi4_x(lam0,mu0);
  phi4(2) = int_phi4_y(lam0,mu0);
  phi4(3) = int_phi4_xy(lam0,mu0);
  // std::cout << "phi4: " << phi4 << '\n';
  phi4_init = phi4;
  std::cout << "phi4_init: " << phi4_init << '\n';
  // Initial values for ODEs
  phi4(0) = kv::complex<intval>(0,0);
  phi4(1) = kv::complex<intval>(0,0);
  phi4(2) = kv::complex<intval>(0,0);
  phi4(3) = kv::complex<intval>(1,0);
  //
  /////// Compute the analytic continuation on the path gamma_{6,1}
  phi4_vec.resize(2*n);
  //
  for (int i = 0; i < n; i++) {
    phi4_vec(i) = phi4(i).real();
    phi4_vec(i+n) = phi4(i).imag();
  }
  // std::cout << "phi_vec: " << phi4_vec << '\n';
  vec_odefunc_segment< kv::complex<intval> > f61(x0,y0,c6x,c6y);
  start4 = 0;
  end4 = 1;
  std::cout << "Computing phi4 along the path gamma_{6,1}..." << '\n';
  r = kv::odelong_maffine(
    f61,
    phi4_vec,
    start4,
    end4,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi4_path61: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi4_path61: Solution incomplitely calculated until t = " << end4 << ".\n";
  } else if (r == 3) {
    std::cout << "phi4_path61: Solution calculated until t = " << end4 << ".\n";
  } else {
    std::cout << "phi4_path61: Solution completely calculated until t = " << end4 << ".\n";
  }
  //
  phi4_path61.resize(n);
  for (int i = 0; i < n; i++) {
    phi4_path61(i).real() = phi4_vec(i);
    phi4_path61(i).imag() = phi4_vec(i+n);
  }
  std::cout << "phi4_path61: " << phi4_path61 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,2}
  vec_odefunc_anticlockwise< kv::complex<intval> > f62(x6,y6,r6);
  start4 = -1;
  end4 = 1;
  std::cout << "Computing phi4 along the path gamma_{6,2}..." << '\n';
  r = kv::odelong_maffine(
    f62,
    phi4_vec,
    start4,
    end4,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi4_path62: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi4_path62: Solution incomplitely calculated until t = " << end4 << ".\n";
  } else if (r == 3) {
    std::cout << "phi4_path62: Solution calculated until t = " << end4 << ".\n";
  } else {
    std::cout << "phi4_path62: Solution completely calculated until t = " << end4 << ".\n";
  }
  //
  phi4_path62.resize(n);
  for (int i = 0; i < n; i++) {
    phi4_path62(i).real() = phi4_vec(i);
    phi4_path62(i).imag() = phi4_vec(i+n);
  }
  std::cout << "phi4_path62: " << phi4_path62 << '\n';
  /////// Compute the analytic continuation on the path gamma_{6,3}
  vec_odefunc_segment< kv::complex<intval> > f63(c6x,c6y,x0,y0);
  start4 = 0;
  end4 = 1;
  std::cout << "Computing phi4 along the path gamma_{6,3}..." << '\n';
  r = kv::odelong_maffine(
    f63,
    phi4_vec,
    start4,
    end4,
    kv::ode_param<dd>().set_verbose(0).set_order(order_of_taylor).set_restart_max(10).set_epsilon((dd)pow(2,-68))
  );
  if (r == 0) {
    std::cout << "phi4_path63: Cannot calculate solution.\n";
  } else if (r == 1) {
    std::cout << "phi4_path63: Solution incomplitely calculated until t = " << end4 << ".\n";
  } else if (r == 3) {
    std::cout << "phi4_path63: Solution calculated until t = " << end4 << ".\n";
  } else {
    std::cout << "phi4_path63: Solution completely calculated until t = " << end4 << ".\n";
  }
  //
  phi4_path63.resize(n);
  for (int i = 0; i < n; i++) {
    phi4_path63(i).real() = phi4_vec(i);
    phi4_path63(i).imag() = phi4_vec(i+n);
  }
  std::cout << "phi4_final: " << phi4_path63 << '\n';
  std::cout << "thread = " << omp_get_thread_num() << '\n';
}
}
//////////////////////////////////////////////////////////////// parallel end ///////////////////////////////////////////////////////
// start to compute monodoromy matrix
// M = verifylss(A,B);
  ub::matrix< kv::complex<intval> > A, B, M;
  ub::matrix<intval> A_r, B_r, M_r;
  A.resize(n,n);
  B.resize(n,n);
  M.resize(n,n);
  A_r.resize(2*n,2*n);
  B_r.resize(2*n,n);
  M_r.resize(2*n,n);

  for (int i  = 0; i < n; i++) {
    A(i,0) = phi1_init(i);
    A(i,1) = phi2_init(i);
    A(i,2) = phi3_init(i);
    A(i,3) = phi4_init(i);
    B(i,0) = phi1_path63(i);
    B(i,1) = phi2_path63(i);
    B(i,2) = phi3_path63(i);
    B(i,3) = phi4_path63(i);
  }
  // A  = Ar + Aci;
  // A_r = [Ar, -Ac; Ac, Ar]
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A_r(i,j) = A(i,j).real();
      A_r(i+N,j) = A(i,j).imag();
      A_r(i,j+N) = -A(i,j).imag();
      A_r(i+N,j+N) = A(i,j).real();
    }
  }
  // B = Br + Bci
  // B_r = [Br;Bc]
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      B_r(i,j) = B(i,j).real();
      B_r(i+N,j) = B(i,j).imag();
    }
  }
  kv::vleq(A_r,B_r,M_r); // Solve M = A\B (aka A_r \ B_r)
  for (int i = 0; i < n; i++) {//kv::complex<intval>(y_vec(i),y_vec(i+N));
    for (int j = 0; j < n; j++) {
      M(i,j) = kv::complex<intval>(M_r(i,j),M_r(i+N,j));
    }
  }
  std::cout << "The monodromy matrix: " << prod(M,A) << "\n";
  return 0;
}
