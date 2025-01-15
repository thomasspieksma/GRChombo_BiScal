//////////////////////////////////////////////
#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT Iteration_max,
                                                              CCTK_INT NmaxLeaver,
                                                              CCTK_COMPLEX init_omega,
                                                              CCTK_REAL mass,
                                                              CCTK_REAL MBH,
                                                              CCTK_REAL ll,
                                                              CCTK_REAL d_omega,
                                                              CCTK_REAL threshold_convergence)
{
   CCTK_INFO("START Leaver method");

   CCTK_COMPLEX omega0,omega1;
   omega0 = init_omega;
   omega1 = init_omega;

   /*
   CCTK_BYTE flag11,flag12,flag21,flag22;
   flag11 = false;
   flag12 = false;
   flag21 = false;
   flag22 = false;
   */
   for(int nn=1;nn<=Iteration_max;++nn)
   {
      CCTK_COMPLEX small_Re,small_Im;
      small_Re = CCTK_Cmplx(d_omega*CCTK_CmplxReal(omega0),0.);
      small_Im = CCTK_Cmplx(0.,d_omega*CCTK_CmplxImag(omega0));

      CCTK_COMPLEX omega_p_small_Re,omega_p_small_Im;
      omega_p_small_Re = CCTK_CmplxAdd(omega0,small_Re);
      omega_p_small_Im = CCTK_CmplxAdd(omega0,small_Im);

      CCTK_COMPLEX res,res_p_small_Re,res_p_small_Im;
      res            = get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(NmaxLeaver,
                                                                     omega0,
                                                                     mass,
                                                                     MBH,
                                                                     ll);
      res_p_small_Re = get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(NmaxLeaver,
                                                                     omega_p_small_Re,
                                                                     mass,
                                                                     MBH,
                                                                     ll);
      res_p_small_Im = get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(NmaxLeaver,
                                                                     omega_p_small_Im,
                                                                     mass,
                                                                     MBH,
                                                                     ll);
      CCTK_REAL d_Re_res_Re,d_Re_res_Im,d_Im_res_Re,d_Im_res_Im;
      d_Re_res_Re = (CCTK_CmplxReal(res_p_small_Re) - CCTK_CmplxReal(res))/CCTK_CmplxReal(small_Re);
      d_Re_res_Im = (CCTK_CmplxImag(res_p_small_Re) - CCTK_CmplxImag(res))/CCTK_CmplxReal(small_Re);
      d_Im_res_Re = (CCTK_CmplxReal(res_p_small_Im) - CCTK_CmplxReal(res))/CCTK_CmplxImag(small_Im);
      d_Im_res_Im = (CCTK_CmplxImag(res_p_small_Im) - CCTK_CmplxImag(res))/CCTK_CmplxImag(small_Im);

      CCTK_REAL Jacobi[2][2];
      Jacobi[0][0] = d_Re_res_Re;
      Jacobi[0][1] = d_Im_res_Re;
      Jacobi[1][0] = d_Re_res_Im;
      Jacobi[1][1] = d_Im_res_Im;
      CCTK_REAL detJ = Jacobi[0][0]*Jacobi[1][1] - Jacobi[1][0]*Jacobi[0][1];
      CCTK_REAL IJacobi[2][2];
      IJacobi[0][0] =  Jacobi[1][1]/detJ;
      IJacobi[0][1] = -Jacobi[0][1]/detJ;
      IJacobi[1][0] = -Jacobi[1][0]/detJ;
      IJacobi[1][1] =  Jacobi[0][0]/detJ;

      omega1 = CCTK_Cmplx(CCTK_CmplxReal(omega0) - (IJacobi[0][0]*CCTK_CmplxReal(res) + IJacobi[0][1]*CCTK_CmplxImag(res)),
                          CCTK_CmplxImag(omega0) - (IJacobi[1][0]*CCTK_CmplxReal(res) + IJacobi[1][1]*CCTK_CmplxImag(res)));
      omega0 = omega1;
   }

   {
      CCTK_COMPLEX res;
      res            = get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(NmaxLeaver,
                                                                     omega0,
                                                                     mass,
                                                                     MBH,
                                                                     ll);

      CCTK_INFO("Leaver method end:");
      CCTK_VInfo(CCTK_THORNSTRING, "omega.Re=>%12.7e  omega.Im=>%12.7e",
                (double)CCTK_CmplxReal(omega0), (double)CCTK_CmplxImag(omega0));
      CCTK_VInfo(CCTK_THORNSTRING, "res.Re=>%12.7e  res.Im=>%12.7e",
                (double)CCTK_CmplxReal(res), (double)CCTK_CmplxImag(res));

      if(threshold_convergence < CCTK_CmplxAbs(res))
      {
         CCTK_WARN(CCTK_WARN_ABORT,"Leaver method did not converge.");
      }
   }
   return omega0;
}

CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                           CCTK_COMPLEX omega,
                                                           CCTK_REAL mass,
                                                           CCTK_REAL MBH,
                                                           CCTK_REAL ll)
{
   CCTK_COMPLEX va1;
   va1 = get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(-1,omega,mass,MBH,ll);

   CCTK_COMPLEX va2;
   va2 = get_Leaver_con_frac_NBI_ScalarKerr_LeaverSchwarzschild(NmaxLeaver,omega,mass,MBH,ll);

   return CCTK_CmplxSub(va1,va2);
}

CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverSchwarzschild(CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll)
{
   CCTK_COMPLEX va = CCTK_Cmplx(1.0,0.0);
   for(int nn=NmaxLeaver;nn>=0;--nn)
   {
      CCTK_COMPLEX alpha = get_alpha_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn-1,omega,mass,MBH,ll);
      CCTK_COMPLEX beta  = get_beta_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn,omega,mass,MBH,ll);
      CCTK_COMPLEX gamma = get_gamma_NBI_ScalarKerr_coef_LeaverSchwarzschild(nn,omega,mass,MBH,ll);

      CCTK_COMPLEX z_num,z_den;
      z_num = CCTK_CmplxMul(alpha,gamma);
      z_den = CCTK_CmplxSub(beta,va);
      va = CCTK_CmplxDiv(z_num,z_den);
   }
   return va;
}

