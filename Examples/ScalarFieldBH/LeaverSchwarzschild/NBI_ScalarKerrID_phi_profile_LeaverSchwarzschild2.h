///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERSCHWARZSCHILD2
#define INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERSCHWARZSCHILD2
//////////////////////////////////////////////

#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_REAL get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(CCTK_REAL x1,
                                                      CCTK_REAL y1,
                                                      CCTK_REAL z1,
                                                      CCTK_REAL tt,
                                                      CCTK_INT NmaxLeaver,
                                                      CCTK_COMPLEX omega,
                                                      CCTK_REAL ini_phase,
                                                      CCTK_REAL mass,
                                                      CCTK_REAL MBH,
                                                      CCTK_INT ll,
                                                      CCTK_INT mm,
                                                      CCTK_COMPLEX *Coef);


CCTK_REAL get_dx_phi_NBI_ScalarKerrID_LeaverSchwazschild2(CCTK_REAL x1,
                                                        CCTK_REAL y1,
                                                        CCTK_REAL z1,
                                                        CCTK_REAL tt,
                                                        CCTK_INT NmaxLeaver,
                                                        CCTK_COMPLEX omega,
                                                        CCTK_REAL ini_phase,
                                                        CCTK_REAL mass,
                                                        CCTK_REAL MBH,
                                                        CCTK_INT ll,
                                                        CCTK_INT mm,
                                                        CCTK_COMPLEX *Coef,
                                                        CCTK_REAL dx);

CCTK_REAL get_dy_phi_NBI_ScalarKerrID_LeaverSchwazschild2(CCTK_REAL x1,
                                                        CCTK_REAL y1,
                                                        CCTK_REAL z1,
                                                        CCTK_REAL tt,
                                                        CCTK_INT NmaxLeaver,
                                                        CCTK_COMPLEX omega,
                                                        CCTK_REAL ini_phase,
                                                        CCTK_REAL mass,
                                                        CCTK_REAL MBH,
                                                        CCTK_INT ll,
                                                        CCTK_INT mm,
                                                        CCTK_COMPLEX *Coef,
                                                        CCTK_REAL dy);

CCTK_REAL get_dz_phi_NBI_ScalarKerrID_LeaverSchwazschild2(CCTK_REAL x1,
                                                        CCTK_REAL y1,
                                                        CCTK_REAL z1,
                                                        CCTK_REAL tt,
                                                        CCTK_INT NmaxLeaver,
                                                        CCTK_COMPLEX omega,
                                                        CCTK_REAL ini_phase,
                                                        CCTK_REAL mass,
                                                        CCTK_REAL MBH,
                                                        CCTK_INT ll,
                                                        CCTK_INT mm,
                                                        CCTK_COMPLEX *Coef,
                                                        CCTK_REAL dz);

CCTK_REAL get_dt_phi_NBI_ScalarKerrID_LeaverSchwazschild2(CCTK_REAL x1,
                                                        CCTK_REAL y1,
                                                        CCTK_REAL z1,
                                                        CCTK_REAL tt,
                                                        CCTK_INT NmaxLeaver,
                                                        CCTK_COMPLEX omega,
                                                        CCTK_REAL ini_phase,
                                                        CCTK_REAL mass,
                                                        CCTK_REAL MBH,
                                                        CCTK_INT ll,
                                                        CCTK_INT mm,
                                                        CCTK_COMPLEX *Coef,
                                                        CCTK_REAL dt);



#endif
