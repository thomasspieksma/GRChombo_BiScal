///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERKERR
#define INCLUDE_HEAD_NBI_SCALARKERRID_PHI_PROFILE_LEAVERKERR
//////////////////////////////////////////////

#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_REAL get_phi_NBI_ScalarKerrID_LeaverKerr(const CCTK_REAL x1,
                                          const CCTK_REAL y1,
                                          const CCTK_REAL z1,
                                          const CCTK_REAL tt,
                                          const CCTK_INT LeaverOrderRadial,
                                          const CCTK_COMPLEX* CoefRadial,
                                          const CCTK_INT LeaverOrderAng,
                                          const CCTK_COMPLEX* CoefAng,
                                          const CCTK_COMPLEX omega,
                                          const CCTK_REAL mass,
                                          const CCTK_REAL MBH,
                                          const CCTK_REAL aBH,
                                          const CCTK_COMPLEX lambda,
                                          const CCTK_INT mm,
                                          const CCTK_REAL rmax_Windowfunc,
                                          const CCTK_REAL rmin_Windowfunc);

CCTK_REAL get_dx_phi_NBI_ScalarKerrID_LeaverKerr(const CCTK_REAL x1,
                                             const CCTK_REAL y1,
                                             const CCTK_REAL z1,
                                             const CCTK_REAL tt,
                                             const CCTK_INT LeaverOrderRadial,
                                             const CCTK_COMPLEX* CoefRadial,
                                             const CCTK_INT LeaverOrderAng,
                                             const CCTK_COMPLEX* CoefAng,
                                             const CCTK_COMPLEX omega,
                                             const CCTK_REAL mass,
                                             const CCTK_REAL MBH,
                                             const CCTK_REAL aBH,
                                             const CCTK_COMPLEX lambda,
                                             const CCTK_INT mm,
                                             const CCTK_REAL rmax_Windowfunc,
                                             const CCTK_REAL rmin_Windowfunc,
                                             const CCTK_REAL dx);

CCTK_REAL get_dy_phi_NBI_ScalarKerrID_LeaverKerr(const CCTK_REAL x1,
                                             const CCTK_REAL y1,
                                             const CCTK_REAL z1,
                                             const CCTK_REAL tt,
                                             const CCTK_INT LeaverOrderRadial,
                                             const CCTK_COMPLEX* CoefRadial,
                                             const CCTK_INT LeaverOrderAng,
                                             const CCTK_COMPLEX* CoefAng,
                                             const CCTK_COMPLEX omega,
                                             const CCTK_REAL mass,
                                             const CCTK_REAL MBH,
                                             const CCTK_REAL aBH,
                                             const CCTK_COMPLEX lambda,
                                             const CCTK_INT mm,
                                             const CCTK_REAL rmax_Windowfunc,
                                             const CCTK_REAL rmin_Windowfunc,
                                             const CCTK_REAL dy);

CCTK_REAL get_dz_phi_NBI_ScalarKerrID_LeaverKerr(const CCTK_REAL x1,
                                             const CCTK_REAL y1,
                                             const CCTK_REAL z1,
                                             const CCTK_REAL tt,
                                             const CCTK_INT LeaverOrderRadial,
                                             const CCTK_COMPLEX* CoefRadial,
                                             const CCTK_INT LeaverOrderAng,
                                             const CCTK_COMPLEX* CoefAng,
                                             const CCTK_COMPLEX omega,
                                             const CCTK_REAL mass,
                                             const CCTK_REAL MBH,
                                             const CCTK_REAL aBH,
                                             const CCTK_COMPLEX lambda,
                                             const CCTK_INT mm,
                                             const CCTK_REAL rmax_Windowfunc,
                                             const CCTK_REAL rmin_Windowfunc,
                                             const CCTK_REAL dz);

CCTK_REAL get_dt_phi_NBI_ScalarKerrID_LeaverKerr(const CCTK_REAL x1,
                                             const CCTK_REAL y1,
                                             const CCTK_REAL z1,
                                             const CCTK_REAL tt,
                                             const CCTK_INT LeaverOrderRadial,
                                             const CCTK_COMPLEX* CoefRadial,
                                             const CCTK_INT LeaverOrderAng,
                                             const CCTK_COMPLEX* CoefAng,
                                             const CCTK_COMPLEX omega,
                                             const CCTK_REAL mass,
                                             const CCTK_REAL MBH,
                                             const CCTK_REAL aBH,
                                             const CCTK_COMPLEX lambda,
                                             const CCTK_INT mm,
                                             const CCTK_REAL rmax_Windowfunc,
                                             const CCTK_REAL rmin_Windowfunc,
                                             const CCTK_REAL dt);

CCTK_REAL get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                            const CCTK_REAL yw,
                                                            const CCTK_REAL zw,
                                                            const CCTK_REAL aBH);
CCTK_REAL get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                           const CCTK_REAL yw,
                                                           const CCTK_REAL zw,
                                                           const CCTK_REAL aBH);
CCTK_REAL get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                         const CCTK_REAL yw,
                                                         const CCTK_REAL zw,
                                                         const CCTK_REAL aBH);

#endif
