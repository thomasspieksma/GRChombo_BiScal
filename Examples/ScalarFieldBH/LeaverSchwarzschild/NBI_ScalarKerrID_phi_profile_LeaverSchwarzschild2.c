#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_spherical_harmonics.h"
#include "NBI_ScalarKerrID_phi_profile_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_phi_profile_LeaverSchwarzschild2.h"

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
                                                       CCTK_COMPLEX *Coef)
{
   const CCTK_COMPLEX imagNum = CCTK_Cmplx(0.,1.);
   const int spin = 0;

   CCTK_REAL rw = get_NBI_ScalarKerrID_radial_coordinate_LeaverSchwazschild(x1,y1,z1);
   if(rw < 1e-5*MBH) { rw = 1e-5*MBH; }
   CCTK_REAL theta = get_NBI_ScalarKerrID_theta_coordinate_LeaverSchwazschild(x1,y1,z1);
   CCTK_REAL varphi = get_NBI_ScalarKerrID_phi_coordinate_LeaverSchwazschild(x1,y1,z1);


   CCTK_COMPLEX rphi =  get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild2(rw,
                                                                               NmaxLeaver,
                                                                               omega,
                                                                               mass,
                                                                               MBH,
                                                                               ll,
                                                                               Coef);


   const CCTK_COMPLEX cmplx1 = CCTK_Cmplx(1.0,0.0);
   const CCTK_COMPLEX cmplxm1= CCTK_Cmplx(-1.0,0.0);
   const CCTK_COMPLEX cmplx2 = CCTK_Cmplx(2.0,0.0);
   const CCTK_COMPLEX omega2 = CCTK_CmplxMul(omega,omega);
   const CCTK_COMPLEX mass2  = CCTK_Cmplx(mass*mass,0.);
   const CCTK_COMPLEX mass2_omega2 = CCTK_CmplxSub(mass2,omega2);
   const CCTK_COMPLEX sqrt_mass2_omega2 = CCTK_CmplxSqrt(mass2_omega2);
   const CCTK_COMPLEX m_sqrt_mass2_omega2 = CCTK_CmplxMul(sqrt_mass2_omega2,cmplxm1);
   const CCTK_COMPLEX mass2_2omega2 = CCTK_CmplxSub(mass2,CCTK_CmplxMul(omega2,cmplx2));
   const CCTK_COMPLEX chi = CCTK_CmplxMul(CCTK_Cmplx(MBH,0),CCTK_CmplxDiv(mass2_2omega2,m_sqrt_mass2_omega2));
   const CCTK_COMPLEX m2Momega = CCTK_Cmplx(-2*MBH*CCTK_CmplxReal(omega),-2*MBH*CCTK_CmplxImag(omega));
   const CCTK_COMPLEX mI2Momega = CCTK_Cmplx(-CCTK_CmplxImag(m2Momega),CCTK_CmplxReal(m2Momega));
   const CCTK_COMPLEX mI2Momega_chi_1 = CCTK_CmplxAdd(mI2Momega,CCTK_CmplxSub(chi,cmplx1));
   const CCTK_COMPLEX logr = CCTK_Cmplx(log(rw),0.);
   const CCTK_COMPLEX r_mI2Momega_chi_1 = CCTK_CmplxExp(CCTK_CmplxMul(mI2Momega_chi_1,logr));
   const CCTK_COMPLEX radial_phi = CCTK_CmplxMul(r_mI2Momega_chi_1,rphi);

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
                                                        CCTK_REAL dx)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1 + 2.*dx,
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
                                                        Coef);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1 + dx,
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
                                                        Coef);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1 - dx,
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
                                                        Coef);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1 - 2.*dx,
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
                                                        Coef);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dx);
}

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
                                                        CCTK_REAL dy)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dy);
}

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
                                                         CCTK_REAL dz)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dz);
}

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
                                                        CCTK_REAL dt)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverSchwazschild2(x1,
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
                                                        Coef);
   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dt);
}

