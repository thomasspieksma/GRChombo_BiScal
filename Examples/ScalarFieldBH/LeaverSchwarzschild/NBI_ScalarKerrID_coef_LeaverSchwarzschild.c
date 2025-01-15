
#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


CCTK_COMPLEX get_alpha_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                               CCTK_COMPLEX omega,
                                                               CCTK_REAL mass,
                                                               CCTK_REAL MBH,
                                                               CCTK_REAL ll)
{
   CCTK_REAL Reomega = CCTK_CmplxReal(omega);
   CCTK_REAL Imomega = CCTK_CmplxImag(omega);

   CCTK_REAL mass2 = pow(mass,2.);

   /*
   CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_Cmplx(mass2 - pow(Reomega,2) + pow(Imomega,2),
                                               -2*Imomega*Reomega));
   */
   //CCTK_REAL Reqq = CCTK_CmplxReal(qq);
   //CCTK_REAL Imqq = CCTK_CmplxImag(qq);

   CCTK_REAL Realpha,Imalpha;

   Realpha = -((2 + nn)*(-4*pow(Imomega,3)*MBH + 
       pow(Imomega,2)*(2 + nn) - 
       4*Imomega*MBH*
        (mass2 - 3*pow(Reomega,2)) + 
       (2 + nn)*(mass2 - pow(Reomega,2))));

   Imalpha = 2*(2 + nn)*Reomega*(-6*pow(Imomega,2)*MBH + 
     Imomega*(2 + nn) + 
     2*MBH*(-mass2 + pow(Reomega,2)));

   return CCTK_Cmplx(Realpha,Imalpha);
}

CCTK_COMPLEX get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                              CCTK_COMPLEX omega,
                                                              CCTK_REAL mass,
                                                              CCTK_REAL MBH,
                                                              CCTK_REAL ll)
{
   CCTK_REAL Reomega = CCTK_CmplxReal(omega);
   CCTK_REAL Imomega = CCTK_CmplxImag(omega);

   CCTK_REAL mass2 = pow(mass,2.);

   CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_Cmplx(mass2 - pow(Reomega,2) + pow(Imomega,2),
                                               -2*Imomega*Reomega));
   CCTK_REAL Reqq = CCTK_CmplxReal(qq);
   CCTK_REAL Imqq = CCTK_CmplxImag(qq);

   CCTK_REAL Rebeta,Imbeta;

   Rebeta = 16*pow(Imomega,4)*pow(MBH,2) + 
   4*pow(mass2,2)*pow(MBH,2) - 
   4*Imomega*mass2*MBH*(3 + 2*nn + 3*MBH*Reqq) - 
   4*pow(Imomega,3)*MBH*
    (3 + 2*nn + 4*MBH*Reqq) + 
   pow(Imomega,2)*
    (5 + ll + pow(ll,2) + 
      20*mass2*pow(MBH,2) + 6*nn + 
      2*pow(nn,2) - 
      48*Imqq*pow(MBH,2)*Reomega - 
      96*pow(MBH,2)*pow(Reomega,2) + 
      12*MBH*Reqq + 8*MBH*nn*Reqq) + 
   mass2*(5 + ll + pow(ll,2) + 2*pow(nn,2) - 
      12*Imqq*pow(MBH,2)*Reomega - 
      20*pow(MBH,2)*pow(Reomega,2) + 
      9*MBH*Reqq + 6*nn*(1 + MBH*Reqq)) + 
   pow(Reomega,2)*
    (-5 - ll - pow(ll,2) - 2*pow(nn,2) + 
      16*Imqq*pow(MBH,2)*Reomega + 
      16*pow(MBH,2)*pow(Reomega,2) - 
      12*MBH*Reqq - 2*nn*(3 + 4*MBH*Reqq)) + 
   4*Imomega*MBH*Reomega*
    (Imqq*(6 + 4*nn) + 
      3*Reomega*(3 + 2*nn + 4*MBH*Reqq));

   Imbeta = -16*pow(Imomega,3)*pow(MBH,2)*
    (Imqq + 4*Reomega) + 
   4*pow(Imomega,2)*MBH*
    (Imqq*(3 + 2*nn) + 
      3*Reomega*(3 + 2*nn + 4*MBH*Reqq)) + 
   MBH*(Imqq*(3 + 2*nn)*
       (3*mass2 - 4*pow(Reomega,2)) + 
      4*mass2*Reomega*(3 + 2*nn + 3*MBH*Reqq) - 
      4*pow(Reomega,3)*(3 + 2*nn + 4*MBH*Reqq))
     - 2*Imomega*(6*Imqq*pow(MBH,2)*
       (mass2 - 4*pow(Reomega,2)) + 
      Reomega*(5 + ll + pow(ll,2) + 
         20*mass2*pow(MBH,2) + 6*nn + 
         2*pow(nn,2) - 
         32*pow(MBH,2)*pow(Reomega,2) + 
         12*MBH*Reqq + 8*MBH*nn*Reqq));
   
  return CCTK_Cmplx(Rebeta,Imbeta);
}

CCTK_COMPLEX get_gamma_NBI_ScalarKerr_coef_LeaverSchwarzschild(CCTK_INT nn,
                                                               CCTK_COMPLEX omega,
                                                               CCTK_REAL mass,
                                                               CCTK_REAL MBH,
                                                               CCTK_REAL ll)
{
   CCTK_REAL Reomega = CCTK_CmplxReal(omega);
   CCTK_REAL Imomega = CCTK_CmplxImag(omega);

   CCTK_REAL mass2 = pow(mass,2.);

   CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_Cmplx(mass2 - pow(Reomega,2) + pow(Imomega,2),
                                               -2*Imomega*Reomega));
   CCTK_REAL Reqq = CCTK_CmplxReal(qq);
   CCTK_REAL Imqq = CCTK_CmplxImag(qq);

   CCTK_REAL Regamma,Imgamma;
 
   Regamma = -8*pow(Imomega,4)*pow(MBH,2) - 
   pow(mass2,2)*pow(MBH,2) + 
   4*pow(Imomega,3)*MBH*
    (1 + nn + 2*MBH*Reqq) - 
   mass2*(1 + 2*nn + pow(nn,2) - 
      4*Imqq*pow(MBH,2)*Reomega - 
      8*pow(MBH,2)*pow(Reomega,2) + 
      2*MBH*Reqq + 2*MBH*nn*Reqq) - 
   pow(Imomega,2)*
    (1 + 8*mass2*pow(MBH,2) + 2*nn + 
      pow(nn,2) - 
      24*Imqq*pow(MBH,2)*Reomega - 
      48*pow(MBH,2)*pow(Reomega,2) + 
      4*MBH*Reqq + 4*MBH*nn*Reqq) + 
   pow(Reomega,2)*
    (1 + 2*nn + pow(nn,2) - 
      8*Imqq*pow(MBH,2)*Reomega - 
      8*pow(MBH,2)*pow(Reomega,2) + 
      4*MBH*Reqq + 4*MBH*nn*Reqq) + 
   4*Imomega*MBH*(mass2*(1 + nn + MBH*Reqq) - 
      Reomega*(2*Imqq*(1 + nn) + 
         3*Reomega*(1 + nn + 2*MBH*Reqq)));

   Imgamma = 8*pow(Imomega,3)*pow(MBH,2)*
    (Imqq + 4*Reomega) - 
   4*pow(Imomega,2)*MBH*
    (Imqq*(1 + nn) + 
      3*Reomega*(1 + nn + 2*MBH*Reqq)) - 
   2*MBH*(Imqq*(1 + nn)*
       (mass2 - 2*pow(Reomega,2)) + 
      2*mass2*Reomega*(1 + nn + MBH*Reqq) - 
      2*pow(Reomega,3)*(1 + nn + 2*MBH*Reqq))\
    + 2*Imomega*(2*Imqq*pow(MBH,2)*
       (mass2 - 6*pow(Reomega,2)) + 
      Reomega*(1 + 8*mass2*pow(MBH,2) + 2*nn + 
         pow(nn,2) - 
         16*pow(MBH,2)*pow(Reomega,2) + 
         4*MBH*Reqq + 4*MBH*nn*Reqq));

   return CCTK_Cmplx(Regamma,Imgamma);
}

