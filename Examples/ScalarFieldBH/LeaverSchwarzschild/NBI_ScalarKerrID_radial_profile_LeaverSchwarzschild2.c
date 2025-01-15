#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild2.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX get_max_NBI_ScalarKerr_LeaverSchwarzschild2(CCTK_INT NmaxLeaver,
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
      phil = get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild2(rw,
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


