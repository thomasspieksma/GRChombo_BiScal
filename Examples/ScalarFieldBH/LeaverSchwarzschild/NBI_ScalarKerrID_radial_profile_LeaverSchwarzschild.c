#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Normalize_coef_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                       CCTK_COMPLEX *Coef,
                                                       CCTK_REAL fac)
{
   for(CCTK_INT nn=0;nn<NmaxLeaver;++nn)
   {
      Coef[nn] = Coef[nn]*fac;
   }
   return;
}

CCTK_COMPLEX get_max_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                     CCTK_COMPLEX omega,
                                                     CCTK_REAL mass,
                                                     CCTK_REAL MBH,
                                                     CCTK_REAL ll,
                                                     CCTK_COMPLEX *Coef,
                                                     CCTK_REAL rmax_search_max,
                                                     CCTK_REAL rmin_search_max,
                                                     CCTK_INT N_grid_search_max)
{
   CCTK_REAL dr = (rmax_search_max - rmin_search_max)/(N_grid_search_max-1);
   CCTK_COMPLEX phi_max = CCTK_Cmplx(0.0, 0.0);
   for(int nn=1;nn<=N_grid_search_max;++nn)
   {
      CCTK_REAL rw = nn*dr + rmin_search_max;
      CCTK_COMPLEX phil;
      phil = get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild(rw,
                                                                   NmaxLeaver,
                                                                   omega,
                                                                   mass,
                                                                   MBH,
                                                                   ll,
                                                                   Coef);

      if(fabs(CCTK_CmplxReal(phil)) > fabs(CCTK_CmplxReal(phi_max)))
      {
         phi_max = phil;
      }
   }

   return phi_max;
}


void get_coef_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_COMPLEX *Coef,
                                                 CCTK_INT NmaxLeaver,
                                                 CCTK_COMPLEX omega,
                                                 CCTK_REAL mass,
                                                 CCTK_REAL MBH,
                                                 CCTK_REAL ll)
{
   Coef[0] = CCTK_Cmplx(1.0,0.0);
   {
      CCTK_COMPLEX alpha = get_alpha_NBI_ScalarKerr_coef_LeaverSchwarzschild(-1,omega,mass,MBH,ll);
      CCTK_COMPLEX beta  = get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(-1,omega,mass,MBH,ll);
      CCTK_COMPLEX beta_alpha = CCTK_CmplxDiv(beta,alpha);
      CCTK_COMPLEX m_beta_alpha = CCTK_Cmplx(-CCTK_CmplxReal(beta_alpha),
                                             -CCTK_CmplxImag(beta_alpha));
      Coef[1] = CCTK_CmplxMul(m_beta_alpha,Coef[0]);
   }

   for(int nn=0;nn<NmaxLeaver-2;++nn)
   {
      CCTK_COMPLEX alpha = get_alpha_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn,omega,mass,MBH,ll);
      CCTK_COMPLEX beta  = get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn,omega,mass,MBH,ll);
      CCTK_COMPLEX gamma = get_gamma_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn,omega,mass,MBH,ll);

      CCTK_COMPLEX beta_alpha = CCTK_CmplxDiv(beta,alpha);
      CCTK_COMPLEX gamma_alpha = CCTK_CmplxDiv(gamma,alpha);

      CCTK_COMPLEX m_beta_alpha = CCTK_Cmplx(-CCTK_CmplxReal(beta_alpha),
                                             -CCTK_CmplxImag(beta_alpha));
      CCTK_COMPLEX m_gamma_alpha = CCTK_Cmplx(-CCTK_CmplxReal(gamma_alpha),
                                              -CCTK_CmplxImag(gamma_alpha));

      Coef[nn+2] = CCTK_CmplxAdd(CCTK_CmplxMul(m_beta_alpha,Coef[nn+1]),
                                 CCTK_CmplxMul(m_gamma_alpha,Coef[nn]));

   }
   return;
}

CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_REAL rw,
                                                                CCTK_INT NmaxLeaver,
                                                                CCTK_COMPLEX omega,
                                                                CCTK_REAL mass,
                                                                CCTK_REAL MBH,
                                                                CCTK_REAL ll,
                                                                CCTK_COMPLEX *Coef)
{
   ////////////
   const CCTK_COMPLEX mass2_omega2 = CCTK_CmplxSub(CCTK_Cmplx(pow(mass,2.),0.0),
                                                   CCTK_CmplxMul(omega,omega));
   const CCTK_COMPLEX sqrt_mass2_omega2 = CCTK_CmplxSqrt(mass2_omega2);
   const CCTK_COMPLEX sqrt_mass2_omega2_r = CCTK_CmplxMul(sqrt_mass2_omega2,
                                                          CCTK_Cmplx(rw,0.0));
   const CCTK_COMPLEX m_sqrt_mass2_omega2_r = CCTK_CmplxMul(sqrt_mass2_omega2_r,
                                                            CCTK_Cmplx(-1.0,0.0));
   const CCTK_COMPLEX exp_m_sqrt_mass2_omega2_r = CCTK_CmplxExp(m_sqrt_mass2_omega2_r);
   ////////////
   const CCTK_REAL r_2M_1 = (rw/(2.*MBH)) - 1.0;
   const CCTK_REAL log_r_2M_1 = log(r_2M_1);
   const CCTK_COMPLEX I2Momega = CCTK_CmplxMul(CCTK_Cmplx(0.0,2.*MBH),
                                               omega);
   const CCTK_COMPLEX r_2M_1_I2Momega = CCTK_CmplxExp(CCTK_CmplxMul(I2Momega,
                                                                    CCTK_Cmplx(log_r_2M_1,0.0)));
  ////////////
  // const CCTK_REAL first = rw - 2.0*MBH;
  // const CCTK_REAL log_first = log(first);
  // const CCTK_COMPLEX factorA = CCTK_CmplxMul(CCTK_Cmplx(0.0,2.*MBH),omega);
  // const CCTK_COMPLEX first_I = CCTK_CmplxExp(CCTK_CmplxMul(factorA,
  //                                                       CCTK_Cmplx(log_first,0.0)));

   ///////////
   const CCTK_REAL mass2 = pow(mass,2.);
   const CCTK_COMPLEX omega2 = CCTK_CmplxMul(omega,omega);
   const CCTK_COMPLEX mass2_2omega2 = CCTK_CmplxSub(CCTK_Cmplx(mass2,0.0),
                                                    CCTK_Cmplx(2.*CCTK_CmplxReal(omega2),
                                                               2.*CCTK_CmplxImag(omega2)));
   const CCTK_COMPLEX mass2_2omega2_ovr_sqrt_mass2_omega2 = CCTK_CmplxDiv(mass2_2omega2,sqrt_mass2_omega2);
   const CCTK_COMPLEX eta = CCTK_CmplxSub(CCTK_Cmplx(-MBH*CCTK_CmplxReal(mass2_2omega2_ovr_sqrt_mass2_omega2),
                                                     -MBH*CCTK_CmplxImag(mass2_2omega2_ovr_sqrt_mass2_omega2)),
                                          I2Momega);
   const CCTK_COMPLEX r_eta = CCTK_CmplxExp(CCTK_CmplxMul(eta,
                                                          CCTK_Cmplx(log(rw),0.0)));
   /////////////
   const CCTK_REAL r_2M_ovr_r = (rw - 2.*MBH)/rw;
   /////////////
   CCTK_COMPLEX sum = CCTK_Cmplx(0.0,0.0);
   for(int nn=0;nn<NmaxLeaver;++nn)
   {
      const CCTK_REAL r_2M_ovr_r_n = pow(r_2M_ovr_r,nn);
      const CCTK_COMPLEX Coef_r_2M_ovr_r_n = CCTK_CmplxMul(Coef[nn],
                                                           CCTK_Cmplx(r_2M_ovr_r_n,0.0));
      sum = CCTK_CmplxAdd(sum,Coef_r_2M_ovr_r_n);
   }

   CCTK_COMPLEX rphi;
   rphi = CCTK_CmplxMul(exp_m_sqrt_mass2_omega2_r,
                        r_2M_1_I2Momega);
   rphi = CCTK_CmplxMul(rphi,r_eta);
   rphi = CCTK_CmplxMul(rphi,sum);
   CCTK_COMPLEX over_rw  = CCTK_Cmplx(1./rw, 0.0);
   //rphi = CCTK_CmplxMul(rphi,first_I);
   CCTK_COMPLEX phi_l = CCTK_CmplxMul(rphi, over_rw);
   return phi_l;
}

CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild2(CCTK_REAL rw,
                                                                    CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll,
                                                                    CCTK_COMPLEX *Coef)
{
   ////////////
   const CCTK_COMPLEX mass2_omega2 = CCTK_CmplxSub(CCTK_Cmplx(pow(mass,2.),0.0),
                                                   CCTK_CmplxMul(omega,omega));
   const CCTK_COMPLEX sqrt_mass2_omega2 = CCTK_CmplxSqrt(mass2_omega2);
   const CCTK_COMPLEX sqrt_mass2_omega2_r = CCTK_CmplxMul(sqrt_mass2_omega2,
                                                          CCTK_Cmplx(rw,0.0));
   const CCTK_COMPLEX m_sqrt_mass2_omega2_r = CCTK_CmplxMul(sqrt_mass2_omega2_r,
                                                            CCTK_Cmplx(-1.0,0.0));
   const CCTK_COMPLEX exp_m_sqrt_mass2_omega2_r = CCTK_CmplxExp(m_sqrt_mass2_omega2_r);
   ////////////
   const CCTK_REAL r_2M_ovr_r = (rw - 2.*MBH)/rw;
   /////////////
   CCTK_COMPLEX sum = CCTK_Cmplx(0.0,0.0);
   for(int nn=0;nn<NmaxLeaver;++nn)
   {
      const CCTK_REAL r_2M_ovr_r_n = pow(r_2M_ovr_r,nn);
      const CCTK_COMPLEX Coef_r_2M_ovr_r_n = CCTK_CmplxMul(Coef[nn],
                                                           CCTK_Cmplx(r_2M_ovr_r_n,0.0));
      sum = CCTK_CmplxAdd(sum,Coef_r_2M_ovr_r_n);
   }

   const CCTK_COMPLEX rphi = CCTK_CmplxMul(sum,exp_m_sqrt_mass2_omega2_r);
   return rphi;
}


