///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_RADIAL_PROFILE_LEAVERKERR2
#define INCLUDE_HEAD_NBI_SCALARKERRID_RADIAL_PROFILE_LEAVERKERR2
//////////////////////////////////////////////

#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


CCTK_COMPLEX get_alpha_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                      const CCTK_COMPLEX omega,
                                                      const CCTK_REAL mass,
                                                      const CCTK_REAL MBH,
                                                      const CCTK_REAL aBH,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm);
CCTK_COMPLEX get_beta_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                     const CCTK_COMPLEX omega,
                                                     const CCTK_REAL mass,
                                                     const CCTK_REAL MBH,
                                                     const CCTK_REAL aBH,
                                                     const CCTK_COMPLEX lambda,
                                                     const CCTK_INT mm);
CCTK_COMPLEX get_gamma_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                      const CCTK_COMPLEX omega,
                                                      const CCTK_REAL mass,
                                                      const CCTK_REAL MBH,
                                                      const CCTK_REAL aBH,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm);
void GetCoef_NBI_ScalarKerr_LeaverKerr(const CCTK_INT LeaverOrderRadial,
                                       const CCTK_COMPLEX lambda,
                                       const CCTK_INT mm,
                                       const CCTK_COMPLEX omega,
                                       const CCTK_REAL mass,
                                       const CCTK_REAL MBH,
                                       const CCTK_REAL aBH,
                                       CCTK_COMPLEX *Coef);

CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverKerr(const CCTK_REAL rw,
                                                          const CCTK_INT LeaverOrderRadial,
                                                          const CCTK_COMPLEX omega,
                                                          const CCTK_REAL mass,
                                                          const CCTK_REAL MBH,
                                                          const CCTK_REAL aBH,
                                                          const CCTK_COMPLEX lambda,
                                                          const CCTK_INT mm,
                                                          const CCTK_COMPLEX *Coef);


CCTK_COMPLEX get_max_NBI_ScalarKerr_LeaverKerr(const CCTK_INT Nrmax,
                                               const CCTK_REAL rmin,
                                               const CCTK_REAL rmax,
                                               const CCTK_INT LeaverOrderRadial,
                                               const CCTK_COMPLEX omega,
                                               const CCTK_REAL mass,
                                               const CCTK_REAL MBH,
                                               const CCTK_REAL aBH,
                                               const CCTK_COMPLEX lambda,
                                               const CCTK_INT mm,
                                               const CCTK_COMPLEX *Coef);
#endif
