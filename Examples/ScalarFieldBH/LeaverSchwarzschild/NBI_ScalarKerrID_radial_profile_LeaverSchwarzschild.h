///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_RADIAL_PROFILE_LEAVERSCHWARZSCHILD
#define INCLUDE_HEAD_NBI_SCALARKERRID_RADIAL_PROFILE_LEAVERSCHWARZSCHILD
//////////////////////////////////////////////
#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Normalize_coef_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                        CCTK_COMPLEX *Coef,
                                                        CCTK_REAL fac);

CCTK_COMPLEX get_max_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                     CCTK_COMPLEX omega,
                                                     CCTK_REAL mass,
                                                     CCTK_REAL MBH,
                                                     CCTK_REAL ll,
                                                     CCTK_COMPLEX *Coef,
                                                     CCTK_REAL rmax_search_max,
                                                     CCTK_REAL rmin_search_max,
                                                     CCTK_INT N_grid_search_max);

void get_coef_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_COMPLEX *Coef,
                                                 CCTK_INT NmaxLeaver,
                                                 CCTK_COMPLEX omega,
                                                 CCTK_REAL mass,
                                                 CCTK_REAL MBH,
                                                 CCTK_REAL ll);

CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_REAL rw,
                                                                CCTK_INT NmaxLeaver,
                                                                CCTK_COMPLEX omega,
                                                                CCTK_REAL mass,
                                                                CCTK_REAL MBH,
                                                                CCTK_REAL ll,
                                                                CCTK_COMPLEX *Coef);


CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild2(CCTK_REAL rw,
                                                                    CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll,
                                                                    CCTK_COMPLEX *Coef);

#endif
