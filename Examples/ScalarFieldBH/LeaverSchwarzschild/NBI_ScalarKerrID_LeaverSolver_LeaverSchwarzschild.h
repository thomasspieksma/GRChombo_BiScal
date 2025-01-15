///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_LEAVERSOLVER_LEAVERSCHWARZSCHILD
#define INCLUDE_HEAD_NBI_SCALARKERRID_LEAVERSOLVER_LEAVERSCHWARZSCHILD
//////////////////////////////////////////////
#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT Iteration_max,
                                                              CCTK_INT NmaxLeaver,
                                                              CCTK_COMPLEX init_omega,
                                                              CCTK_REAL mass,
                                                              CCTK_REAL MBH,
                                                              CCTK_REAL ll,
                                                              CCTK_REAL d_omega,
                                                              CCTK_REAL threshold_convergence);

CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                           CCTK_COMPLEX omega,
                                                           CCTK_REAL mass,
                                                           CCTK_REAL MBH,
                                                           CCTK_REAL ll);

CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll);

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

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT Iteration_max,
                                                              CCTK_INT NmaxLeaver,
                                                              CCTK_COMPLEX init_omega,
                                                              CCTK_REAL mass,
                                                              CCTK_REAL MBH,
                                                              CCTK_REAL ll,
                                                              CCTK_REAL d_omega,
                                                              CCTK_REAL threshold_convergence);


CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                           CCTK_COMPLEX omega,
                                                           CCTK_REAL mass,
                                                           CCTK_REAL MBH,
                                                           CCTK_REAL ll);


CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll);

#endif
