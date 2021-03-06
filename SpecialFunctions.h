// [06/2018] Misha Mikhasenko, mikhail.mikhasenko@gmail.com

#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>
#ifndef SPECFUNC_H
#define SPECFUNC_H

namespace SpecialFunc {
//------------------------------------------------------------------------------//

typedef unsigned int uint;

const static double logfact[24] = {0.0, 0.0,
                                     6.93147180559945309e-1, 1.79175946922805500e00,
                                     3.17805383034794562e00, 4.78749174278204599e00,
                                     6.57925121201010100e00, 8.52516136106541430e00,
                                     1.06046029027452502e01, 1.28018274800814696e01,
                                     1.51044125730755153e01, 1.75023078458738858e01,
                                     1.99872144956618861e01, 2.25521638531234229e01,
                                     2.51912211827386815e01, 2.78992713838408916e01,
                                     3.06718601060806728e01, 3.35050734501368889e01,
                                     3.63954452080330536e01, 3.93398841871994940e01,
                                     4.23356164607534850e01, 4.53801388984769080e01,
                                     4.84711813518352239e01, 5.16066755677643736e01};

inline int abs(int x) { return (x > 0) ? x : - x; }
inline double ClebschGordan(uint j1, int m1, uint j2, int m2, uint j, int m) {
        return ((j1-j2+m) % 4 == 0 ? 1.0 : -1.0) * sqrt(j+1) *
               gsl_sf_coupling_3j(j1, j2, j, m1, m2, -m);
}

// The implemented expression can be find in Wikipedia
//     https://en.wikipedia.org/wiki/Jacobi_polynomials
inline double JacobiPols(uint n, uint a, uint b, double x) {
        // here comes the calculation of log(n!)

        if (n+a >= 24 || n+b >= 24) {
                std::cerr << "Error: j is too high, please check the implementation of jacobi polynomials!\n";
                return 0.0;
        }

        double ls = log((1.0-x)/2.0);
        double lc = log((1.0+x)/2.0);
        double res = 0.0;
        for (uint s = 0; s <= n; s++) {
                double logs = logfact[n+a] + logfact[n+b]-logfact[n-s]-logfact[a+s]-logfact[s]-logfact[n+b-s];
                double args = s*ls + (n-s)*lc;
                res += (s % 2 == 0 ? 1.0 : -1.0) * exp(logs+args);
        }
        return res;
}

// The reference  for the relation between WignerD and Jacobi polynomials is
// Eq. (3.74) of L. Biedenharn, J. Louck, and P. Carruthers, Angular Momentum in Quantum Physics: Theory and Application
// see also (B1) of the FormalismII paper.
// j, m1, m2 are doubled
inline double WignerDhat(uint j, int m1, int m2, double z) {
        double factor = ((abs(m1-m2)+m1-m2)) % 8 == 0 ? 1.0 : -1.0;
        int am1 = abs(m1), am2 = abs(m2);
        if (j % 2 != am1 % 2 || j % 2 != am2 % 2)  { std::cerr << "Error: j, m1, m2 not compatible\n"; return 0.; }
        int M = (am1 > am2) ? am1 : am2;
        int N = (am1 < am2) ? am1 : am2;
        return factor/pow(2,.5*M)*
               sqrt(gsl_sf_gamma((j-M)/2+1.)*gsl_sf_gamma((j+M)/2+1.)/(gsl_sf_gamma((j-N)/2+1.)*gsl_sf_gamma((j+N)/2+1.)))*
               JacobiPols((j-M)/2, abs(m1-m2)/2,abs(m1+m2)/2, z);
}

inline double WignerDhat(uint j, int m1, int m2, bool eta, double z) {
        int am1 = abs(m1), am2 = abs(m2);
        int M = (am1 > am2) ? am1 : am2;
        double factor = (m2 - M) % 4 == 0 ? 1.0 : -1.0;
        if (!eta) factor = -factor;
        // std::cout << WignerDhat(j, m1, m2, z) << std::endl;
        return WignerDhat(j, m1, m2, z) + factor*WignerDhat(j, -m1, m2, z);

}

// j, m1, m2 are doubled
inline double WignerD(uint j, int m1, int m2, double z) {
        double hat = WignerDhat(j, m1, m2, z);
        double xi = pow(sqrt(1-z),abs(m1-m2)/2)*pow(sqrt((1+z)),abs(m1+m2)/2);
        // std::cout << "hat " << hat << " xi " << xi << std::endl;
        return hat*xi;
}


}

#endif
