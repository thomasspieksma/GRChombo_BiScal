///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERSCHWARZSCHILD
#define INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERSCHWARZSCHILD
//////////////////////////////////////////////

#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_REAL get_phi_NBI_ScalarKerrID_LeaverSchwazschild(CCTK_REAL x1,
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
                                                      CCTK_REAL rmax_Windowfunc,
                                                      CCTK_REAL rmin_Windowfunc);


CCTK_REAL get_dx_phi_NBI_ScalarKerrID_LeaverSchwazschild(CCTK_REAL x1,
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
                                                        CCTK_REAL rmax_Windowfunc,
                                                        CCTK_REAL rmin_Windowfunc,
                                                        CCTK_REAL dx);

CCTK_REAL get_dy_phi_NBI_ScalarKerrID_LeaverSchwazschild(CCTK_REAL x1,
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
                                                        CCTK_REAL rmax_Windowfunc,
                                                        CCTK_REAL rmin_Windowfunc,
                                                        CCTK_REAL dy);

CCTK_REAL get_dz_phi_NBI_ScalarKerrID_LeaverSchwazschild(CCTK_REAL x1,
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
                                                        CCTK_REAL rmax_Windowfunc,
                                                        CCTK_REAL rmin_Windowfunc,
                                                        CCTK_REAL dz);

CCTK_REAL get_dt_phi_NBI_ScalarKerrID_LeaverSchwazschild(CCTK_REAL x1,
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
                                                        CCTK_REAL rmax_Windowfunc,
                                                        CCTK_REAL rmin_Windowfunc,
                                                        CCTK_REAL dt);




CCTK_REAL get_NBI_ScalarKerrID_radial_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw);
CCTK_REAL get_NBI_ScalarKerrID_theta_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw);
CCTK_REAL get_NBI_ScalarKerrID_phi_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw);

#endif
