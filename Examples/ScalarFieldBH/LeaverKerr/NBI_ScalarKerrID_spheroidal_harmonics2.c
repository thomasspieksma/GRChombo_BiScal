#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "NBI_ScalarKerrID_spheroidal_harmonics2.h"

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_alpha(const CCTK_INT nn,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm,
                                                      const CCTK_COMPLEX gg) {
  const CCTK_REAL re = -2*(2 + nn)*(2 + nn + fabs(mm));
  const CCTK_REAL im = 0.0;
  return CCTK_Cmplx(re,im);
}

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_beta(const CCTK_INT nn,
                                                     const CCTK_COMPLEX lambda,
                                                     const CCTK_INT mm,
                                                     const CCTK_COMPLEX gg) {

  const CCTK_REAL Relambda = CCTK_CmplxReal(lambda);
  const CCTK_REAL Imlambda = CCTK_CmplxImag(lambda);
  const CCTK_REAL Regg     = CCTK_CmplxReal(gg);
  const CCTK_REAL Imgg     = CCTK_CmplxImag(gg);

  const CCTK_REAL re = 2 + pow(mm,2) + 3*nn + pow(nn,2) + Imgg*(6 + 4*nn) - Relambda +
                      (3 + 2*Imgg + 2*nn)*fabs(mm);
  const CCTK_REAL im = -Imlambda - 2*(3 + 2*nn)*Regg - 2*Regg*fabs(mm);

  return CCTK_Cmplx(re,im);
}

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_gamma(const CCTK_INT nn,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm,
                                                      const CCTK_COMPLEX gg) {
  const CCTK_REAL Relambda = CCTK_CmplxReal(lambda);
  const CCTK_REAL Imlambda = CCTK_CmplxImag(lambda);
  const CCTK_REAL Regg     = CCTK_CmplxReal(gg);
  const CCTK_REAL Imgg     = CCTK_CmplxImag(gg);

  const CCTK_REAL re = -2*Imgg*(1 + nn + fabs(mm));
  const CCTK_REAL im = 2*Regg*(1 + nn + fabs(mm));

  return CCTK_Cmplx(re,im);
}

void nbi_scalarkerr_SpheroidalHarmonics_GetCoef(const CCTK_INT LeaverOrderAng,
                                                const CCTK_INT lambda,
                                                const CCTK_INT mm,
                                                const CCTK_COMPLEX gg,
                                                CCTK_COMPLEX *Coef) {
  Coef[0] = 1;

  {
    const CCTK_COMPLEX alpha = nbi_scalarkerr_SpheroidalHarmonics_alpha(-1,lambda,mm,gg);
    const CCTK_COMPLEX beta  = nbi_scalarkerr_SpheroidalHarmonics_beta(-1,lambda,mm,gg);
    const CCTK_COMPLEX va1   = CCTK_CmplxDiv(beta,alpha);
    const CCTK_COMPLEX va2   = CCTK_CmplxMul(va1,Coef[0]);
    Coef[1] = CCTK_Cmplx(-CCTK_CmplxReal(va2),-CCTK_CmplxImag(va2));
  }

  const CCTK_COMPLEX m1 = CCTK_Cmplx(-1.,0.);

  for(int nn=0;nn<LeaverOrderAng-2;++nn){
    const CCTK_COMPLEX alpha = nbi_scalarkerr_SpheroidalHarmonics_alpha(nn,lambda,mm,gg);
    const CCTK_COMPLEX beta  = nbi_scalarkerr_SpheroidalHarmonics_beta(nn,lambda,mm,gg);
    const CCTK_COMPLEX gamma = nbi_scalarkerr_SpheroidalHarmonics_gamma(nn,lambda,mm,gg);

    const CCTK_COMPLEX beta_alpha = CCTK_CmplxDiv(beta,alpha);
    const CCTK_COMPLEX gamma_alpha = CCTK_CmplxDiv(gamma,alpha);

    const CCTK_COMPLEX v1 = CCTK_CmplxMul(beta_alpha,Coef[nn+1]);
    const CCTK_COMPLEX v2 = CCTK_CmplxMul(gamma_alpha,Coef[nn]);

    Coef[nn+2] = CCTK_CmplxMul(m1,CCTK_CmplxAdd(v1,v2));

      //cout << "alpha,beta,gamma " << alpha.Re << "," << alpha.Im << "," << beta.Re << "," << beta.Im <<"," <<  gamma.Re << "," << gamma.Im << endl;
   }

  const CCTK_INT imax = 100;
  const CCTK_REAL dth = M_PI/imax;
  CCTK_REAL sum=0.0;
  for(int ii=0;ii<imax;++ii) {
    const CCTK_REAL th = (ii + 0.5)*dth;
    const CCTK_REAL xw = cos(th);

    const CCTK_COMPLEX exp_I_gamma
      = CCTK_CmplxExp(CCTK_CmplxMul(CCTK_Cmplx(0.0,xw),gg));
    const CCTK_REAL one_m_x = 1.0 - xw;
    const CCTK_REAL one_p_x = 1.0 + xw;
    const CCTK_REAL one_m_x_absm_2 = pow(one_m_x,fabs(mm)*0.5);
    const CCTK_REAL one_p_x_absm_2 = pow(one_p_x,fabs(mm)*0.5);

    CCTK_COMPLEX fac = CCTK_Cmplx(0.0,0.0);
    for(int nn=0;nn<LeaverOrderAng;++nn) {
      fac = CCTK_CmplxAdd(fac,CCTK_CmplxMul(Coef[nn],pow(one_p_x,nn)));
    }

    const CCTK_COMPLEX S = CCTK_CmplxMul(fac,CCTK_CmplxMul(exp_I_gamma,one_m_x_absm_2*one_p_x_absm_2));
    sum += pow(CCTK_CmplxAbs(S),2)*sin(th);
  }

  const CCTK_REAL integral = sum*dth;
  const CCTK_REAL sqrt_int = sqrt(integral);
  for(int nw=0;nw<LeaverOrderAng;++nw) {
    Coef[nw] = CCTK_CmplxDiv(Coef[nw],sqrt_int);
  }
  return;
}

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_GetValue(const CCTK_INT LeaverOrderAng,
                                                         const CCTK_INT lambda,
                                                         const CCTK_INT mm,
                                                         const CCTK_COMPLEX gg,
                                                         CCTK_COMPLEX *Coef,
                                                         const CCTK_REAL th,
                                                         const CCTK_REAL ph) {
  const CCTK_REAL xw = cos(th);

  const CCTK_COMPLEX exp_I_gamma
      = CCTK_CmplxExp(CCTK_CmplxMul(CCTK_Cmplx(0.0,xw),gg));
  const CCTK_REAL one_m_x = 1.0 - xw;
  const CCTK_REAL one_p_x = 1.0 + xw;
  const CCTK_REAL one_m_x_absm_2 = pow(one_m_x,fabs(mm)*0.5);
  const CCTK_REAL one_p_x_absm_2 = pow(one_p_x,fabs(mm)*0.5);

  CCTK_COMPLEX fac = CCTK_Cmplx(0.0,0.0);
  for(int nn=0;nn<LeaverOrderAng;++nn) {
    fac = CCTK_CmplxAdd(fac,CCTK_CmplxMul(Coef[nn],pow(one_p_x,nn)));
  }

  return CCTK_CmplxMul(CCTK_CmplxMul(fac,
                                      CCTK_CmplxMul(exp_I_gamma,one_m_x_absm_2*one_p_x_absm_2))
                      ,CCTK_CmplxExp(CCTK_Cmplx(0.0,mm*ph)));
}

