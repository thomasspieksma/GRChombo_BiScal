#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_spherical_harmonics.h"
#include "NBI_ScalarKerrID_phi_profile_LeaverSchwarzschild.h"

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
                                                      CCTK_REAL rmin_Windowfunc)
{
   const CCTK_COMPLEX imagNum = CCTK_Cmplx(0.,1.);
   const int spin = 0;

   CCTK_REAL rw = get_NBI_ScalarKerrID_radial_coordinate_LeaverSchwazschild(x1,y1,z1);
   CCTK_REAL theta = get_NBI_ScalarKerrID_theta_coordinate_LeaverSchwazschild(x1,y1,z1);
   CCTK_REAL varphi = get_NBI_ScalarKerrID_phi_coordinate_LeaverSchwazschild(x1,y1,z1);


   CCTK_COMPLEX radial_phi =  get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild(rw,
                                                                                   NmaxLeaver,
                                                                                   omega,
                                                                                   mass,
                                                                                   MBH,
                                                                                   ll,
                                                                                   Coef);

  // CCTK_COMPLEX radial_phi = CCTK_Cmplx(radial_phi_re, 0.0);
   CCTK_COMPLEX first = CCTK_Cmplx(rw - 2.0*MBH,0.0);
   CCTK_COMPLEX log_first = CCTK_CmplxLog(first);
  // CCTK_COMPLEX log_first_c = CCTK_Cmplx(log_first, 0.0);
   CCTK_COMPLEX factorA = CCTK_CmplxMul(CCTK_Cmplx(2.*MBH, 0.0),omega);
   factorA = CCTK_CmplxMul(factorA, CCTK_Cmplx(-1., 0.0));
   factorA = CCTK_CmplxMul(factorA, imagNum);
   CCTK_COMPLEX additionalI = CCTK_CmplxExp(CCTK_CmplxMul(factorA,log_first));
   //printf("r: %5.2e \n", rw);
   //printf("Factor: %5.2e \n", additionalI);
   //printf(%5.2e\n,additional);
   radial_phi = CCTK_CmplxMul(radial_phi,additionalI);
  
   CCTK_REAL lambda;
   lambda = 1.;
   if(rmin_Windowfunc < rw && rw < rmax_Windowfunc) {
      lambda = 0.50 * ( 1 + tanh( 1.0/(rmin_Windowfunc - rw) - 1.0/(rw - rmax_Windowfunc) ) );
   }
   CCTK_COMPLEX lambdac;
   lambdac = CCTK_Cmplx(lambda,0.0);

   radial_phi = CCTK_CmplxMul(radial_phi,lambdac);

   if(rw <= rmin_Windowfunc) {
      radial_phi = CCTK_Cmplx(0.0,0.0);
   }

   CCTK_REAL reY,imY;
   nbi_scalarkerr_qlm_Multipole_SphericalHarmonic(spin,ll,mm,
                                                  theta,varphi,
                                                  &reY, &imY);

   CCTK_COMPLEX Ylm = CCTK_Cmplx(reY,imY);
 
   CCTK_COMPLEX I_omega_t = CCTK_CmplxMul(imagNum,
                                          CCTK_Cmplx(CCTK_CmplxReal(omega)*tt,
                                                     CCTK_CmplxImag(omega)*tt));
   CCTK_COMPLEX I_omega_t_phase = CCTK_CmplxAdd(I_omega_t,
                                                CCTK_Cmplx(0.0,ini_phase));

   CCTK_COMPLEX exp_I_omega_t_phase = CCTK_CmplxExp(I_omega_t_phase);


   CCTK_COMPLEX phil;
   phil = CCTK_CmplxMul(exp_I_omega_t_phase,
                        CCTK_CmplxMul(Ylm,
                                      radial_phi));
 
   return CCTK_CmplxReal(phil);
}

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
                                                        CCTK_REAL dx)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1 + 2.*dx,
                                                        y1,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1 + dx,
                                                        y1,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1 - dx,
                                                        y1,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1 - 2.*dx,
                                                        y1,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dx);
}

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
                                                        CCTK_REAL dy)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1 + 2.*dy,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1 + dy,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1 - dy,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1 - 2.*dy,
                                                        z1,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dy);
}

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
                                                         CCTK_REAL dz)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1 + 2.*dz,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1 + dz,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1 - dz,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1 - 2.*dz,
                                                        tt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dz);
}

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
                                                        CCTK_REAL dt)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1,
                                                        tt + 2.*dt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1,
                                                        tt + dt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1,
                                                        tt - dt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild(x1,
                                                        y1,
                                                        z1,
                                                        tt - 2.*dt,
                                                        NmaxLeaver,
                                                        omega,
                                                        ini_phase,
                                                        mass,
                                                        MBH,
                                                        ll,
                                                        mm,
                                                        Coef,
                                                        rmax_Windowfunc,
                                                        rmin_Windowfunc);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dt);
}
CCTK_REAL get_NBI_ScalarKerrID_radial_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw)
{   
    return sqrt(xw*xw + yw*yw + zw*zw);
}
CCTK_REAL get_NBI_ScalarKerrID_theta_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw)
{   
    const CCTK_REAL small_val = 1e-10; 
    CCTK_REAL rw        = sqrt(xw*xw + yw*yw + zw*zw);
    CCTK_REAL thetaw    = acos(zw/rw);
    if(fabs(rw) < small_val) {
       thetaw = 0.;
    }
    return thetaw;
}
CCTK_REAL get_NBI_ScalarKerrID_phi_coordinate_LeaverSchwazschild(CCTK_REAL xw,CCTK_REAL yw,CCTK_REAL zw)
{
    const CCTK_REAL small_val = 1e-10;
    //CCTK_REAL rw        = sqrt(xw*xw + yw*yw + zw*zw);
    CCTK_REAL phiw      = atan(yw/xw);
    if(xw < 0) {
        phiw = atan(yw/xw) + M_PI;
    }
    if(fabs(xw) < small_val && yw >= 0.) {
        phiw = M_PI/2.;
    }
    if(fabs(xw) < small_val && yw < 0.) {
        phiw = -M_PI/2.;
    }
    return phiw;
}

