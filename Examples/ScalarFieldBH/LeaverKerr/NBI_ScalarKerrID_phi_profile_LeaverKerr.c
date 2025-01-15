#include <math.h>
#include <stdbool.h>

//#include "NBI_ScalarKerrID_coef_LeaverKerr.h"
//#include "NBI_ScalarKerrID_LeaverSolver_LeaverKerr.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverKerr2.h"
#include "NBI_ScalarKerrID_spheroidal_harmonics2.h"
#include "NBI_ScalarKerrID_phi_profile_LeaverKerr.h"

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
                                              const CCTK_REAL rmin_Windowfunc)
{
   const CCTK_REAL rw = get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(x1,y1,z1,aBH);
   const CCTK_REAL theta = get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(x1,y1,z1,aBH);
   const CCTK_REAL varphi = get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(x1,y1,z1,aBH);
   const CCTK_COMPLEX radial_phi =
     get_radial_profile_NBI_ScalarKerr_LeaverKerr(rw,
                                                  LeaverOrderRadial,
                                                  omega,
                                                  mass,
                                                  MBH,
                                                  aBH,
                                                  lambda,
                                                  mm,
                                                  CoefRadial);

   const CCTK_COMPLEX gg = CCTK_CmplxSqrt(CCTK_CmplxSub(mass*mass,CCTK_CmplxMul(omega,omega)));

   const CCTK_COMPLEX Slm =
     nbi_scalarkerr_SpheroidalHarmonics_GetValue(LeaverOrderAng,
                                                 lambda,
                                                 mm,
                                                 gg,
                                                 CoefAng,
                                                 theta,
                                                 varphi);


  //if(-0.1<x1&&x1<0.1&&-2.6<y1&&y1<-2.4&&0.74<z1&&z1<0.76) {
    //printf("rw,theta,varphi = (%e,%e,%e)\n",rw,theta,varphi);
    //printf("radial_phi = (%e,%e)\n",CCTK_CmplxReal(radial_phi),CCTK_CmplxImag(radial_phi));
    //printf("Slm = (%e,%e)\n",CCTK_CmplxReal(Slm),CCTK_CmplxImag(Slm));
  //}

  const CCTK_REAL rp = MBH + sqrt(MBH*MBH - aBH*aBH);
  const CCTK_REAL rm = MBH - sqrt(MBH*MBH - aBH*aBH);
  const CCTK_COMPLEX AA = CCTK_CmplxMul(CCTK_Cmplx(0.0,2.0),CCTK_CmplxMul(omega,MBH*rp/(rp - rm)));
  const CCTK_COMPLEX BB = CCTK_CmplxMul(CCTK_Cmplx(0.0,-2.0),CCTK_CmplxMul(omega,MBH*rm/(rp-rm)));
  const CCTK_COMPLEX CC = CCTK_CmplxMul(CCTK_Cmplx(0.0,1.0),mm*aBH/(rm - rp));

  const CCTK_COMPLEX fac1 = CCTK_CmplxExp(CCTK_CmplxMul(CCTK_Cmplx(0.0,-tt),omega));
  const CCTK_COMPLEX fac2 = CCTK_CmplxExp(CCTK_CmplxMul(AA,log(fabs(rw - rp))));
  const CCTK_COMPLEX fac3 = CCTK_CmplxExp(CCTK_CmplxMul(BB,log(fabs(rw - rm))));
  const CCTK_COMPLEX fac4 = CCTK_CmplxExp(CCTK_CmplxMul(CC,log(fabs((rw - rp)/(rw - rm)))));

   CCTK_REAL wind;
   wind = 1.;
   if(rmin_Windowfunc < rw && rw < rmax_Windowfunc) {
      wind = 0.50 * ( 1 + tanh( 1.0/(rmin_Windowfunc - rw) - 1.0/(rw - rmax_Windowfunc) ) );
   }
   CCTK_COMPLEX windc;
   windc = CCTK_Cmplx(wind,0.0);

   CCTK_COMPLEX radial_phil = CCTK_CmplxMul(radial_phi,windc);

   if(rw <= rmin_Windowfunc) {
      radial_phil = CCTK_Cmplx(0.0,0.0);
   }

  const CCTK_COMPLEX phil = CCTK_CmplxMul(fac1,
                                          CCTK_CmplxMul(fac2,
                                          CCTK_CmplxMul(fac3,
                                          CCTK_CmplxMul(radial_phil,Slm))));


   return CCTK_CmplxReal(phil);
}

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
                                             const CCTK_REAL dx){
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1 + 2.*dx,
                                                y1,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1 + dx,
                                                y1,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);



   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1 - dx,
                                                y1,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1 - 2.*dx,
                                                y1,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

  // const CCTK_REAL dx_phi_profile = (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dx);
  //if(fabs(dx_phi_profile) > 10.0 && fabs(z1) < 1.0e-1) {
  //  printf("rw = (%e,%e,%e,%e)\n",
  //          get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(x1 + 2.*dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(x1 + dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(x1 - dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(x1 - 2.*dx,y1,z1,aBH));

  //  printf("theta = (%e,%e,%e,%e)\n",
  //          get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(x1 + 2.*dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(x1 + dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(x1 - dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(x1 - 2.*dx,y1,z1,aBH));

  //  printf("varphi = (%e,%e,%e,%e)\n",
  //          get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(x1 + 2.*dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(x1 + dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(x1 - dx,y1,z1,aBH),
  //          get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(x1 - 2.*dx,y1,z1,aBH));

  //  printf("phi = (%e,%e,%e,%e)\n",phi_p2,phi_p1,phi_m1,phi_m2);
  //  printf("(x,y,z) = (%e,%e,%e)\n",x1,y1,z1);
  //}

   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dx);
}

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
                                             const CCTK_REAL dy)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1 + 2.*dy,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1 + dy,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1 - dy,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1 - 2.*dy,
                                                z1,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dy);
}

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
                                             const CCTK_REAL dz)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1,
                                                z1 + 2.*dz,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);


   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1,
                                                z1 + dz,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1,
                                                z1 - dz,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                y1,
                                                z1 - 2.*dz,
                                                tt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dz);
}

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
                                             const CCTK_REAL dt)
{
   CCTK_REAL phi_p2,phi_p1,phi_m1,phi_m2;
   phi_p2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                        y1,
                                                        z1,
                                                        tt + 2.*dt,
                                                LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_p1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                        y1,
                                                        z1,
                                                        tt + dt,
                                              LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);


   phi_m1 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                        y1,
                                                        z1,
                                                        tt - dt,
                                              LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);

   phi_m2 = get_phi_NBI_ScalarKerrID_LeaverKerr(x1,
                                                        y1,
                                                        z1,
                                                        tt - 2.*dt,
                                              LeaverOrderRadial,
                                                CoefRadial,
                                                LeaverOrderAng,
                                                CoefAng,
                                                omega,
                                                mass,
                                                MBH,
                                                aBH,
                                                lambda,
                                                mm,
                                                rmax_Windowfunc,
                                                rmin_Windowfunc);


   return (-phi_p2 + 8.*phi_p1 - 8.*phi_m1 + phi_m2)/(12.*dt);
}
CCTK_REAL get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                            const CCTK_REAL yw,
                                                            const CCTK_REAL zw,
                                                            const CCTK_REAL aBH)
{
    const CCTK_REAL x2_y2_z2 = xw*xw + yw*yw + zw*zw;
    return sqrt(0.5*(x2_y2_z2 - aBH*aBH) + sqrt(pow(x2_y2_z2 - aBH*aBH,2.0) + 4.0*aBH*aBH*zw*zw));
}
CCTK_REAL get_NBI_ScalarKerrID_theta_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                           const CCTK_REAL yw,
                                                           const CCTK_REAL zw,
                                                           const CCTK_REAL aBH)
{
    const CCTK_REAL small_val = 1e-10;
    const CCTK_REAL rw        = get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(xw,yw,zw,aBH);
    CCTK_REAL thetaw    = acos(zw/rw);
    if(fabs(rw) < small_val) {
       thetaw = 0.;
    }
    return thetaw;
}
CCTK_REAL get_NBI_ScalarKerrID_phi_coordinate_LeaverKerr(const CCTK_REAL xw,
                                                         const CCTK_REAL yw,
                                                         const CCTK_REAL zw,
                                                         const CCTK_REAL aBH)
{
    const CCTK_REAL small_val = 1e-10;
    const CCTK_REAL rw        = get_NBI_ScalarKerrID_radial_coordinate_LeaverKerr(xw,yw,zw,aBH);
    const CCTK_REAL XX = rw*xw + aBH*yw;
    const CCTK_REAL YY = rw*yw - aBH*xw;
    //CCTK_REAL phiw      = atan((rw*yw - aBH*xw)/(rw*xw + aBH*yw));
    //if(xw < 0.0) {
    //    phiw = phiw + M_PI;
    //}
    //if(fabs(rw*xw + aBH*yw) < small_val) {
    //    phiw = M_PI/2.;
    //}
    //if(fabs(xw) < small_val && yw < 0. && fabs(aBH) < small_val) {
    //    phiw = -M_PI/2.;
    //}
    const CCTK_REAL phiw = atan2(YY,XX);
    return phiw;
}

