///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_LEAVERSOLVER_LEAVERKERR
#define INCLUDE_HEAD_NBI_SCALARKERRID_LEAVERSOLVER_LEAVERKERR
//////////////////////////////////////////////
#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT Iteration_max,
                                                     const CCTK_INT NmaxLeaver,
                                                     const CCTK_COMPLEX init_omega,
                                                     const CCTK_REAL mass,
                                                     const CCTK_REAL MBH,
                                                     const CCTK_REAL aBH,
                                                     const CCTK_REAL ll,
                                                     const CCTK_REAL mm,
                                                     const CCTK_REAL d_omega,
                                                     const CCTK_REAL threshold_convergence);

CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                  const CCTK_COMPLEX omega,
                                                  const CCTK_REAL mass,
                                                  const CCTK_REAL MBH,
                                                  const CCTK_REAL aBH,
                                                  const CCTK_REAL ll,
                                                  const CCTK_REAL mm);

CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                           const CCTK_COMPLEX omega,
                                                           const CCTK_REAL mass,
                                                           const CCTK_REAL MBH,
                                                           const CCTK_REAL aBH,
                                                           const CCTK_REAL ll,
                                                           const CCTK_REAL mm);

void get_coef_NBI_ScalarKerr_LeaverKerr(CCTK_COMPLEX *Coef,
                                        const CCTK_INT NmaxLeaver,
                                        const CCTK_COMPLEX omega,
                                        const CCTK_REAL mass,
                                        const CCTK_REAL MBH,
                                        const CCTK_REAL aBH,
                                        const CCTK_REAL ll,
                                        const CCTK_REAL mm);

//CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverKerr(const CCTK_REAL rw,
//                                                          const CCTK_INT NmaxLeaver,
//                                                          const CCTK_COMPLEX omega,
//                                                          const CCTK_REAL mass,
//                                                          const CCTK_REAL MBH,
//                                                          const CCTK_REAL aBH,
//                                                          const CCTK_REAL ll,
//                                                          const CCTK_REAL mm,
//                                                          CCTK_COMPLEX *Coef);

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT Iteration_max,
                                                     const CCTK_INT NmaxLeaver,
                                                     const CCTK_COMPLEX init_omega,
                                                     const CCTK_REAL mass,
                                                     const CCTK_REAL MBH,
                                                     const CCTK_REAL aBH,
                                                     const CCTK_REAL ll,
                                                     const CCTK_REAL mm,
                                                     const CCTK_REAL d_omega,
                                                     const CCTK_REAL threshold_convergence);


CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                  const CCTK_COMPLEX omega,
                                                  const CCTK_REAL mass,
                                                  const CCTK_REAL MBH,
                                                  const CCTK_REAL aBH,
                                                  const CCTK_REAL ll,
                                                  const CCTK_REAL mm);


CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                           const CCTK_COMPLEX omega,
                                                           const CCTK_REAL mass,
                                                           const CCTK_REAL MBH,
                                                           const CCTK_REAL aBH,
                                                           const CCTK_REAL ll,
                                                           const CCTK_REAL mm);

#endif
