//////////////////////////////////////////////
#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_LeaverSolver_LeaverKerr.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverKerr2.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX Leaver_solver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT Iteration_max,
                                                     const CCTK_INT NmaxLeaver,
                                                     const CCTK_COMPLEX init_omega,
                                                     const CCTK_REAL mass,
                                                     const CCTK_REAL MBH,
                                                     const CCTK_REAL aBH,
                                                     const CCTK_REAL ll,
                                                     const CCTK_REAL mm,
                                                     const CCTK_REAL d_omega,
                                                     const CCTK_REAL threshold_convergence)
{
   CCTK_INFO("START Leaver method");

   CCTK_COMPLEX omega0,omega1;
   omega0 = init_omega;
   omega1 = init_omega;

   CCTK_REAL d_omega_l = d_omega;

   CCTK_REAL fac = 1e-1;

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
      small_Re = CCTK_Cmplx(d_omega_l*CCTK_CmplxReal(omega0),0.);
      small_Im = CCTK_Cmplx(0.,d_omega_l*CCTK_CmplxImag(omega0));

      CCTK_COMPLEX omega_p_small_Re,omega_p_small_Im;
      omega_p_small_Re = CCTK_CmplxAdd(omega0,small_Re);
      omega_p_small_Im = CCTK_CmplxAdd(omega0,small_Im);

      CCTK_VInfo(CCTK_THORNSTRING, "omega0.Re=>%12.7e  omega0.Im=>%12.7e",
                (double)CCTK_CmplxReal(omega0), (double)CCTK_CmplxImag(omega0));

      CCTK_COMPLEX res,res_p_small_Re,res_p_small_Im;
      res            = get_Leaver_NBI_ScalarKerr_LeaverKerr(NmaxLeaver,
                                                            omega0,
                                                            mass,
                                                            MBH,
                                                            aBH,
                                                            ll,
                                                            mm);
      CCTK_VInfo(CCTK_THORNSTRING, "res.Re=>%12.7e  res.Im=>%12.7e",
                (double)CCTK_CmplxReal(res), (double)CCTK_CmplxImag(res));

      //exit(0);

      res_p_small_Re = get_Leaver_NBI_ScalarKerr_LeaverKerr(NmaxLeaver,
                                                            omega_p_small_Re,
                                                            mass,
                                                            MBH,
                                                            aBH,
                                                            ll,
                                                            mm);
      CCTK_VInfo(CCTK_THORNSTRING, "res_p_small_Re.Re=>%12.7e  res_p_small_Re.Im=>%12.7e",
                (double)CCTK_CmplxReal(res_p_small_Re), (double)CCTK_CmplxImag(res_p_small_Re));

      res_p_small_Im = get_Leaver_NBI_ScalarKerr_LeaverKerr(NmaxLeaver,
                                                            omega_p_small_Im,
                                                            mass,
                                                            MBH,
                                                            aBH,
                                                            ll,
                                                            mm);
      CCTK_VInfo(CCTK_THORNSTRING, "res_p_small_Im.Re=>%12.7e  res_p_small_Im.Im=>%12.7e",
                (double)CCTK_CmplxReal(res_p_small_Im), (double)CCTK_CmplxImag(res_p_small_Im));



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

      omega1 = CCTK_Cmplx(CCTK_CmplxReal(omega0) - fac*(IJacobi[0][0]*CCTK_CmplxReal(res) + IJacobi[0][1]*CCTK_CmplxImag(res)),
                          CCTK_CmplxImag(omega0) - fac*(IJacobi[1][0]*CCTK_CmplxReal(res) + IJacobi[1][1]*CCTK_CmplxImag(res)));

      const CCTK_COMPLEX dev_omega = CCTK_CmplxSub(omega1,omega0);
      omega0 = omega1;
      if(fabs(CCTK_CmplxImag(dev_omega)) < d_omega_l) {
        d_omega_l = CCTK_CmplxAbs(dev_omega);
      }
   }

   {
      CCTK_COMPLEX res;
      res            = get_Leaver_NBI_ScalarKerr_LeaverKerr(NmaxLeaver,
                                                            omega0,
                                                            mass,
                                                            MBH,
                                                            aBH,
                                                            ll,
                                                            mm);

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

CCTK_COMPLEX get_Leaver_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                  const CCTK_COMPLEX omega,
                                                  const CCTK_REAL mass,
                                                  const CCTK_REAL MBH,
                                                  const CCTK_REAL aBH,
                                                  const CCTK_REAL ll,
                                                  const CCTK_REAL mm)
{
  CCTK_COMPLEX va1;
  const CCTK_COMPLEX beta0 = get_beta_NBI_ScalarKerr_coef_LeaverKerr(0,omega,mass,MBH,aBH,ll,mm);
  const CCTK_COMPLEX alph0 = get_alpha_NBI_ScalarKerr_coef_LeaverKerr(0,omega,mass,MBH,aBH,ll,mm);

  //printf("alph0,beta0 = %12.7e, %12.7e, %12.7e, %12.7e\n",CCTK_CmplxReal(alph0),CCTK_CmplxImag(alph0),CCTK_CmplxReal(beta0),CCTK_CmplxImag(beta0));
   va1 = CCTK_CmplxDiv(beta0,alph0);

   CCTK_COMPLEX va2;
   va2 = get_Leaver_con_frac_NBI_ScalarKerr_LeaverKerr(NmaxLeaver,omega,mass,MBH,aBH,ll,mm);

   return CCTK_CmplxSub(va1,va2);
}

CCTK_COMPLEX get_Leaver_con_frac_NBI_ScalarKerr_LeaverKerr(const CCTK_INT NmaxLeaver,
                                                           const CCTK_COMPLEX omega,
                                                           const CCTK_REAL mass,
                                                           const CCTK_REAL MBH,
                                                           const CCTK_REAL aBH,
                                                           const CCTK_REAL ll,
                                                           const CCTK_REAL mm)
{
   CCTK_COMPLEX va = CCTK_Cmplx(-1.0,0.0);
   for(CCTK_INT nn=NmaxLeaver;nn>=1;--nn)
   {
      const CCTK_COMPLEX alpha
        = get_alpha_NBI_ScalarKerr_coef_LeaverKerr(nn,omega,mass,MBH,aBH,ll,mm);
      const CCTK_COMPLEX beta
        = get_beta_NBI_ScalarKerr_coef_LeaverKerr(nn,omega,mass,MBH,aBH,ll,mm);
      const CCTK_COMPLEX gamma
        = get_gamma_NBI_ScalarKerr_coef_LeaverKerr(nn,omega,mass,MBH,aBH,ll,mm);

      const CCTK_COMPLEX alpha_va = CCTK_CmplxMul(alpha,va);

      const CCTK_COMPLEX z_numera = gamma;
      const CCTK_COMPLEX z_denom  = CCTK_CmplxSub(beta,alpha_va);

      va = CCTK_CmplxDiv(z_numera,z_denom);
      //printf("nn=%d,alpha=(%e,%e),beta=(%e,%e),gamma=(%e,%e),va=(%e,%e)\n",nn,CCTK_CmplxReal(alpha),CCTK_CmplxImag(alpha),CCTK_CmplxReal(beta),CCTK_CmplxImag(beta),CCTK_CmplxReal(gamma),CCTK_CmplxImag(gamma),CCTK_CmplxReal(va),CCTK_CmplxImag(va));
   }
   return va;
}

