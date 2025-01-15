///////////////////////////////////////////////↲
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_COEF_LEAVERSCHWARZSCHILD
#define INCLUDE_HEAD_NBI_SCALARKERRID_COEF_LEAVERSCHWARZSCHILD
//////////////////////////////////////////////↲

#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


CCTK_COMPLEX get_alpha_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                               CCTK_COMPLEX omega,
                                                               CCTK_REAL mass,
                                                               CCTK_REAL MBH,
                                                               CCTK_REAL ll);

CCTK_COMPLEX get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                              CCTK_COMPLEX omega,
                                                              CCTK_REAL mass,
                                                              CCTK_REAL MBH,
                                                              CCTK_REAL ll);

CCTK_COMPLEX get_gamma_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                               CCTK_COMPLEX omega,
                                                               CCTK_REAL mass,
                                                               CCTK_REAL MBH,
                                                               CCTK_REAL ll);

#endif
