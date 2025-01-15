///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_SPHERICAL_HARMONICS
#define INCLUDE_HEAD_NBI_SCALARKERRID_SPHERICAL_HARMONICS
//////////////////////////////////////////////


#include <math.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

double nbi_scalarkerr_qlm_factorial(int n);
double nbi_scalarkerr_qlm_combination(int n, int m);
int nbi_scalarkerr_qlm_imin(int a, int b);
int nbi_scalarkerr_qlm_imax(int a, int b);
void nbi_scalarkerr_qlm_Multipole_SphericalHarmonic(int s, int l, int m,
                                                    CCTK_REAL th, CCTK_REAL ph,
                                                    CCTK_REAL *reY, CCTK_REAL *imY);
#endif
