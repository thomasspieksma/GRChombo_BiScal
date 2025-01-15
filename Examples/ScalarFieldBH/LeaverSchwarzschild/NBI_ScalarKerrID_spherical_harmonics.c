#include <math.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "NBI_ScalarKerrID_spherical_harmonics.h"

double nbi_scalarkerr_qlm_factorial(int n)
{
  double returnval = 1;
  for (int i = n; i >= 1; i--)
  {
    returnval *= i;
  }
  return returnval;
}

double nbi_scalarkerr_qlm_combination(int n, int m)
{
  // Binomial coefficient is undefined if these conditions do not hold
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);

  return nbi_scalarkerr_qlm_factorial(n) / (nbi_scalarkerr_qlm_factorial(m) * nbi_scalarkerr_qlm_factorial(n-m));
}

int nbi_scalarkerr_qlm_imin(int a, int b)
{
  return a < b ? a : b;
}

int nbi_scalarkerr_qlm_imax(int a, int b)
{
  return a > b ? a : b;
}

void nbi_scalarkerr_qlm_Multipole_SphericalHarmonic(int s, int l, int m,
                                                    CCTK_REAL th, CCTK_REAL ph,
                                                    CCTK_REAL *reY, CCTK_REAL *imY)
{
  const CCTK_REAL PI = acos(-1.0);
//  assert(s == -2 && l == 2 && m == 2);
//  *reY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * cos(2*ph);
//  *imY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * sin(2*ph);
  double all_coeff = 0, sum = 0;
  all_coeff = pow(-1.0, m);
  all_coeff *= sqrt(nbi_scalarkerr_qlm_factorial(l+m)*nbi_scalarkerr_qlm_factorial(l-m)*(2*l+1) / (4.*PI*nbi_scalarkerr_qlm_factorial(l+s)*nbi_scalarkerr_qlm_factorial(l-s)));
  sum = 0.;
  for (int i = nbi_scalarkerr_qlm_imax(m - s, 0); i <= nbi_scalarkerr_qlm_imin(l + m, l - s); i++)
  {
    double sum_coeff = nbi_scalarkerr_qlm_combination(l-s, i) * nbi_scalarkerr_qlm_combination(l+s, i+s-m);
    sum += sum_coeff * pow(-1.0, l-i-s) * pow(cos(th/2.), 2 * i + s - m) *
      pow(sin(th/2.), 2*(l-i)+m-s);
  }
  *reY = all_coeff*sum*cos(m*ph);
  *imY = all_coeff*sum*sin(m*ph);
}

