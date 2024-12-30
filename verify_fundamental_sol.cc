/*
 * verify_fundamental_sol.cc:
 * This code provides the values of fundametal solutions at the base point (lam, mu) = (2^-10,2^-10)
 *
 * Complile options are as follows (for instance):
 * c++ -I.. -O3 -DNDEBUG verify_fundamental_sol.cc
 * c++ -I.. -O3 -DNDEBUG -DKV_USE_AVX512 -mavx512f -DKV_USE_TPFMA verify_fundamental_sol.cc
 * 
 * written by A. Takayasu	Dec. 28 2024
 */
// 
#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/complex.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
// 
namespace ub = boost::numeric::ublas;
typedef double dd;
// typedef kv::dd dd;
typedef kv::interval<dd> intval;
// 
#include "initial_values.hpp"
// 
const int N = 4; // dimention of ODE (# of complex variables)
const kv::complex<intval> lam0(pow(2, -10)), mu0(pow(2, -10));
// 
int main()
{
    // std::cout << "Max threads: " << omp_get_max_threads() << '\n';
    std::cout.precision(17);
    int n = N; // # of complex variables
    // int r;
    std::cout << "beta_N: " << beta_N << '\n';
    std::cout << "gamm_N: " << gamm_N << '\n';
    std::cout << "eta1_N: " << eta1_N << '\n';
    std::cout << "eta1/2_N: " << eta2_N << '\n';
    std::cout << "theta1_N: " << theta1_N << '\n';
    std::cout << "theta1/2_N: " << theta2_N << '\n';
    std::cout << "iota_N: " << iota_N << '\n';
    std::cout << "nu1_N: " << nu1_N << '\n';
    std::cout << "nu2_N: " << nu2_N << '\n';
    std::cout << "xi1_N: " << xi1_N << '\n';
    std::cout << "xi2_N: " << xi2_N << '\n';
    std::cout << "sig1_N: " << sig1_N << '\n';
    std::cout << "sig2_N: " << sig2_N << '\n';
    std::cout << "C1: " << C1 << '\n';
    std::cout << "delta1: " << delta1 << '\n';
    std::cout << "delta2: " << delta2 << '\n';
    std::cout << "delta3: " << delta3 << '\n';
    std::cout << "delta4: " << delta4 << '\n';
    std::cout << "delta5^1: " << delta51 << '\n';
    std::cout << "delta5^(1/2): " << delta52 << '\n';
    std::cout << "delta6: " << delta6 << '\n';
    std::cout << "delta7: " << delta7 << '\n';
    std::cout << "delta8: " << delta8 << '\n';
    std::cout << "delta9: " << delta9 << '\n';
    std::cout << "delta10^1: " << delta101 << '\n';
    std::cout << "delta10^(1/2): " << delta102 << '\n';
    std::cout << "delta11^1: " << delta111 << '\n';
    std::cout << "delta11^(1/2): " << delta112 << '\n';
    std::cout << "delta12: " << delta12 << '\n';
    std::cout << "delta13: " << delta13 << '\n';
    std::cout << "delta14: " << delta14 << '\n';
    std::cout << "delta15: " << delta15 << '\n';
    std::cout << "delta16^1: " << delta161 << '\n';
    std::cout << "delta16^(1/2): " << delta162 << '\n';
    std::cout << "delta17^1: " << delta171 << '\n';
    std::cout << "delta17^(1/2): " << delta172 << '\n';
    std::cout << "delta18^1: " << delta181 << '\n';
    std::cout << "delta18^(1/2): " << delta182 << '\n';

    /////// Initial value via fundamental solution phi1, phi2, phi3, phi4
    ub::vector<kv::complex<intval>> phi1; //, phi1_init;
    ub::vector<kv::complex<intval>> phi2; //, phi2_init;
    ub::vector<kv::complex<intval>> phi3; //, phi3_init;
    ub::vector<kv::complex<intval>> phi4; //, phi4_init;
    phi1.resize(n);
    phi2.resize(n);
    phi3.resize(n);
    phi4.resize(n);

    phi1(0) = int_phi1(lam0, mu0);
    phi1(1) = int_phi1_x(lam0, mu0);
    phi1(2) = int_phi1_y(lam0, mu0);
    phi1(3) = int_phi1_xy(lam0, mu0);
    std::cout << "phi1: " << phi1 << '\n';

    phi2(0) = int_phi2(lam0, mu0);
    phi2(1) = int_phi2_x(lam0, mu0);
    phi2(2) = int_phi2_y(lam0, mu0);
    phi2(3) = int_phi2_xy(lam0, mu0);
    std::cout << "phi2: " << phi2 << '\n';

    phi3(0) = int_phi3(lam0, mu0);    // this takes more than 5mins
    phi3(1) = int_phi3_x(lam0, mu0);  // this takes more than 5mins
    phi3(2) = int_phi3_y(lam0, mu0);  // this takes more than 5mins
    phi3(3) = int_phi3_xy(lam0, mu0); // this takes more than 5mins
    std::cout << "phi3: " << phi3 << '\n';

    phi4(0) = int_phi4(lam0, mu0);
    phi4(1) = int_phi4_x(lam0, mu0);
    phi4(2) = int_phi4_y(lam0, mu0);
    phi4(3) = int_phi4_xy(lam0, mu0);
    std::cout << "phi4: " << phi4 << '\n';

    return 0;
}
