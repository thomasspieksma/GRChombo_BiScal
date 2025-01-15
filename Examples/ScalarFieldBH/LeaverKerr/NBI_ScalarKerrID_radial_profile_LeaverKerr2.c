#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_radial_profile_LeaverKerr2.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_COMPLEX get_alpha_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                      const CCTK_COMPLEX omega,
                                                      const CCTK_REAL mass,
                                                      const CCTK_REAL MBH,
                                                      const CCTK_REAL aBH,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm) {
  const CCTK_REAL n = (double) nn;

  const CCTK_REAL omegaRe = CCTK_CmplxReal(omega);
  const CCTK_REAL omegaIm = CCTK_CmplxImag(omega);

  const CCTK_REAL b = sqrt(1.0-aBH*aBH);
  const CCTK_REAL a = aBH;
  const CCTK_REAL m = (double) mm;

  const CCTK_REAL c0re = 1.0 + 2.0*(1.0 + 1.0/b)*omegaIm;
  const CCTK_REAL c0im = (a*m - 2.0*(1.0 + b)*omegaRe)/b;

  const CCTK_REAL re = (1.0 + n)*(c0re + n);
  const CCTK_REAL im = (1.0 + n)*c0im;


   /*
   CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_Cmplx(mass2 - pow(Reomega,2) + pow(Imomega,2),
                                               -2*Imomega*Reomega));
   */
   //CCTK_REAL Reqq = CCTK_CmplxReal(qq);
   //CCTK_REAL Imqq = CCTK_CmplxImag(qq);

   return CCTK_Cmplx(re,im);
}

CCTK_COMPLEX get_beta_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                     const CCTK_COMPLEX omega,
                                                     const CCTK_REAL mass,
                                                     const CCTK_REAL MBH,
                                                     const CCTK_REAL aBH,
                                                     const CCTK_COMPLEX lambda,
                                                     const CCTK_INT mm) {

  const CCTK_REAL lambdalmRe = CCTK_CmplxReal(lambda);
  const CCTK_REAL lambdalmIm = CCTK_CmplxImag(lambda);
  const CCTK_REAL omegaRe    = CCTK_CmplxReal(omega);
  const CCTK_REAL omegaIm    = CCTK_CmplxImag(omega);
  const CCTK_REAL a = aBH;
  const CCTK_REAL b = sqrt(1.0-aBH*aBH);
  const CCTK_REAL m = (double) mm;
  const CCTK_REAL n = (double) nn;

  CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_CmplxSub(mass*mass,CCTK_CmplxMul(omega,omega)));
  if(CCTK_CmplxReal(qq) > 0.0) {
    qq = CCTK_CmplxMul(-1.0,qq);
  }
  const CCTK_REAL qRe = CCTK_CmplxReal(qq);
  const CCTK_REAL qIm = CCTK_CmplxImag(qq);
  //printf("n=%d,lambda = %e,%e\n",nn,lambdalmRe,lambdalmIm);

  //std::cout << "lambda = " << lambdalmRe << "," << lambdalmIm << std::endl;

  const CCTK_REAL c1Re = -4.0 - (4.0*omegaIm)/b + 4.0*b*qRe + (2.0*(pow(omegaIm,2.0)*qRe +
        qRe*(-pow(omegaRe,2.0) + pow(qIm,2.0) + pow(qRe,2.0)) -
        2.0*omegaIm*(omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0))))/(pow(qIm,2.0) + pow(qRe,2.0));
  const CCTK_REAL c1Im = (-2.0*a*m + 4.0*omegaRe)/b + 4.0*b*qIm + (2.0*(-(pow(omegaIm,2.0)*qIm) + pow(omegaRe,2.0)*qIm -
        2.0*omegaIm*omegaRe*qRe + 2.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
        qIm*(pow(qIm,2.0) + pow(qRe,2.0))))/(pow(qIm,2.0) + pow(qRe,2.0));

  const CCTK_REAL c3Re = -1.0 + (2.0*pow(b,2.0)*(-pow(omegaIm,2.0) + pow(omegaRe,2.0) + 2.0*omegaRe*qIm + pow(qIm,2.0) + qRe +
         2.0*omegaIm*qRe - pow(qRe,2.0))*(pow(qIm,2.0) + pow(qRe,2.0)) +
      pow(a,2.0)*b*(-pow(qIm,4.0) + pow(qRe,4.0)) +
      a*m*(pow(omegaIm,2.0)*qIm - pow(omegaRe,2.0)*qIm + 2.0*omegaIm*omegaRe*qRe -
         2.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) - (1.0 + 2.0*b)*qIm*(pow(qIm,2.0) + pow(qRe,2.0))) +
      2.0*(pow(omegaIm,3.0)*qRe + omegaIm*(pow(qIm,2.0)*(-1.0 + qRe) +
            qRe*(-3.0*pow(omegaRe,2.0) + (-1.0 + qRe)*qRe)) -
         pow(omegaIm,2.0)*(3.0*omegaRe*qIm + 2.0*(pow(qIm,2.0) + pow(qRe,2.0))) +
         omegaRe*(pow(omegaRe,2.0)*qIm + 2.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
            qIm*(pow(qIm,2.0) + pow(qRe,2.0)))) +
      b*(2.0*pow(omegaRe,3.0)*qIm + 2.0*pow(omegaIm,3.0)*qRe +
         pow(omegaIm,2.0)*(-6.0*omegaRe*qIm - 6.0*pow(qIm,2.0) + qRe - 6.0*pow(qRe,2.0)) +
         6.0*omegaRe*qIm*(pow(qIm,2.0) + pow(qRe,2.0)) +
         (-lambdalmRe + 2.0*pow(qIm,2.0) + qRe - 2.0*pow(qRe,2.0))*(pow(qIm,2.0) + pow(qRe,2.0)) +
         pow(omegaRe,2.0)*(6.0*pow(qIm,2.0) + qRe*(-1.0 + 6.0*qRe)) -
         2.0*omegaIm*(omegaRe*qIm + 3.0*pow(omegaRe,2.0)*qRe - (-1.0 + 3.0*qRe)*(pow(qIm,2.0) + pow(qRe,2.0)))))/
    (b*(pow(qIm,2.0) + pow(qRe,2.0)));

  const CCTK_REAL c3Im = (-2.0*pow(omegaIm,3.0)*qIm + pow(omegaIm,2.0)*(a*m - 6.0*omegaRe)*qRe +
     2.0*pow(b,2.0)*(qIm + 2.0*omegaIm*(omegaRe + qIm) - 2.0*omegaRe*qRe - 2.0*qIm*qRe)*
      (pow(qIm,2.0) + pow(qRe,2.0)) - (a*m - 2.0*omegaRe)*
      (-(pow(qIm,2.0)*(-1.0 + qRe)) + qRe*(pow(omegaRe,2.0) + qRe - pow(qRe,2.0))) +
     2.0*omegaIm*(3.0*pow(omegaRe,2.0)*qIm + 4.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
        qIm*(pow(qIm,2.0) + pow(qRe,2.0)) - a*m*(omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0))) +
     b*(-2.0*pow(omegaIm,3.0)*qIm + pow(omegaRe,2.0)*qIm + 2.0*pow(omegaRe,3.0)*qRe -
        pow(omegaIm,2.0)*(qIm + 6.0*omegaRe*qRe) - 2.0*omegaRe*(-1.0 + 3.0*qRe)*(pow(qIm,2.0) + pow(qRe,2.0)) +
        (-lambdalmIm + qIm + 2.0*a*m*qRe - 4.0*qIm*qRe + 2.0*pow(a,2.0)*qIm*qRe)*(pow(qIm,2.0) + pow(qRe,2.0)) +
        2.0*omegaIm*(3.0*pow(omegaRe,2.0)*qIm + 6.0*omegaRe*pow(qIm,2.0) + omegaRe*qRe*(-1.0 + 6.0*qRe) +
           3.0*qIm*(pow(qIm,2.0) + pow(qRe,2.0)))))/(b*(pow(qIm,2.0) + pow(qRe,2.0)));

  const CCTK_REAL re = c3Re + (2.0 + c1Re - 2.0*n)*n;
  const CCTK_REAL im = c3Im + c1Im*n;

  //printf("n=%d,c1 = %e,%e c3 = %e,%e\n",nn,c1Re,c1Im,c3Re,c3Im);
  //printf("n=%d,re,im = %e,%e\n",nn,re,im);

  return CCTK_Cmplx(re,im);

}

CCTK_COMPLEX get_gamma_NBI_ScalarKerr_coef_LeaverKerr(const CCTK_INT nn,
                                                      const CCTK_COMPLEX omega,
                                                      const CCTK_REAL mass,
                                                      const CCTK_REAL MBH,
                                                      const CCTK_REAL aBH,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm) {

  const CCTK_REAL n = (double) nn;

  const CCTK_REAL omegaRe = CCTK_CmplxReal(omega);
  const CCTK_REAL omegaIm = CCTK_CmplxImag(omega);

  const CCTK_REAL mass2 = pow(mass,2.0);

  const CCTK_REAL b = sqrt(1.0-aBH*aBH);

  const CCTK_REAL q2Re = mass2 - omegaRe*omegaRe + omegaIm*omegaIm;
  const CCTK_REAL q2Im = -2.0*omegaRe*omegaIm;
  const CCTK_REAL sqrt_abs_q2 = sqrt(q2Re*q2Re + q2Im*q2Im);
  CCTK_REAL qRe,qIm;
  if(q2Re > 0.0) {
    qRe = -sqrt(( q2Re + sqrt_abs_q2)*0.5);
    qIm = (q2Im>0) ? -sqrt((-q2Re + sqrt_abs_q2)*0.5)
                   : +sqrt((-q2Re + sqrt_abs_q2)*0.5);
  } else {
    qRe = -sqrt( -q2Im*q2Im*0.5/(q2Re - sqrt_abs_q2));
    qIm = (q2Im>0) ? -sqrt( q2Im*q2Im*0.5/(q2Re + sqrt_abs_q2))
                   : +sqrt( q2Im*q2Im*0.5/(q2Re + sqrt_abs_q2));
  }

  const CCTK_REAL a = aBH;
  const CCTK_REAL m = mm;

  const CCTK_REAL c2Re = 3.0 + (2.0*omegaIm)/b - (2.0*(pow(omegaIm,2.0)*qRe + qRe*(-pow(omegaRe,2.0) + pow(qIm,2.0) + pow(qRe,2.0)) -
        omegaIm*(2.0*omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0))))/(pow(qIm,2.0) + pow(qRe,2.0));
  const CCTK_REAL c2Im = (a*m - 2.0*omegaRe)/b - (2.0*(-(pow(omegaIm,2.0)*qIm) + pow(omegaRe,2.0)*qIm + omegaRe*pow(qIm,2.0) +
        pow(qIm,3.0) - 2.0*omegaIm*omegaRe*qRe + omegaRe*pow(qRe,2.0) + qIm*pow(qRe,2.0)))/
    (pow(qIm,2.0) + pow(qRe,2.0));

  const CCTK_REAL c4Re = (-(b*(pow(omegaIm,4.0)*(pow(qIm,2.0) - pow(qRe,2.0)) + pow(omegaRe,4.0)*(pow(qIm,2.0) - pow(qRe,2.0)) +
          2.0*pow(omegaRe,3.0)*qIm*(pow(qIm,2.0) + pow(qRe,2.0)) +
          2.0*pow(omegaRe,2.0)*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0) +
          2.0*omegaRe*qIm*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0) +
          (pow(qIm,2.0) - pow(qRe,2.0))*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0) +
          2.0*pow(omegaIm,3.0)*qRe*(4.0*omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0)) +
          2.0*omegaIm*qRe*(-4.0*pow(omegaRe,3.0)*qIm - 3.0*pow(omegaRe,2.0)*(pow(qIm,2.0) + pow(qRe,2.0)) +
             pow(pow(qIm,2.0) + pow(qRe,2.0),2.0)) -
          2.0*pow(omegaIm,2.0)*(3.0*pow(omegaRe,2.0)*(pow(qIm,2.0) - pow(qRe,2.0)) +
             3.0*omegaRe*qIm*(pow(qIm,2.0) + pow(qRe,2.0)) + pow(pow(qIm,2.0) + pow(qRe,2.0),2.0)))) +
     (pow(qIm,2.0) + pow(qRe,2.0))*(a*m*(-(pow(omegaIm,2.0)*qIm) + pow(omegaRe,2.0)*qIm -
           2.0*omegaIm*omegaRe*qRe + 2.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
           qIm*(pow(qIm,2.0) + pow(qRe,2.0))) -
        2.0*(pow(omegaIm,3.0)*qRe + omegaIm*qRe*(-3.0*pow(omegaRe,2.0) + pow(qIm,2.0) + pow(qRe,2.0)) -
           pow(omegaIm,2.0)*(3.0*omegaRe*qIm + 2.0*(pow(qIm,2.0) + pow(qRe,2.0))) +
           omegaRe*(pow(omegaRe,2.0)*qIm + 2.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
              qIm*(pow(qIm,2.0) + pow(qRe,2.0))))))/(b*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0));
  const CCTK_REAL c4Im = -((2.0*b*(pow(omegaIm,4.0)*qIm*qRe - 3.0*pow(omegaIm,2.0)*omegaRe*qRe*
           (2.0*omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0)) -
          pow(omegaIm,3.0)*(2.0*omegaRe*pow(qIm,2.0) + pow(qIm,3.0) - 2.0*omegaRe*pow(qRe,2.0) +
             qIm*pow(qRe,2.0)) - qRe*(-(pow(omegaRe,4.0)*qIm) -
             pow(omegaRe,3.0)*(pow(qIm,2.0) + pow(qRe,2.0)) +
             omegaRe*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0) + qIm*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0)) +
          omegaIm*(2.0*pow(omegaRe,3.0)*(pow(qIm,2.0) - pow(qRe,2.0)) +
             3.0*pow(omegaRe,2.0)*qIm*(pow(qIm,2.0) + pow(qRe,2.0)) +
             2.0*omegaRe*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0) + qIm*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0)))\
        + (pow(qIm,2.0) + pow(qRe,2.0))*(-2.0*pow(omegaIm,3.0)*qIm +
          pow(omegaIm,2.0)*(a*m - 6.0*omegaRe)*qRe -
          (a*m - 2.0*omegaRe)*qRe*(pow(omegaRe,2.0) - pow(qIm,2.0) - pow(qRe,2.0)) +
          2.0*omegaIm*(3.0*pow(omegaRe,2.0)*qIm + 4.0*omegaRe*(pow(qIm,2.0) + pow(qRe,2.0)) +
             qIm*(pow(qIm,2.0) + pow(qRe,2.0)) - a*m*(omegaRe*qIm + pow(qIm,2.0) + pow(qRe,2.0)))))/
     (b*pow(pow(qIm,2.0) + pow(qRe,2.0),2.0)));

  const CCTK_REAL re = c4Re + n*(-3.0 + c2Re + n);
  const CCTK_REAL im = c4Im + c2Im*n;

  //printf("n=%d,c2 = %e,%e c4 = %e,%e\n",nn,c2Re,c2Im,c4Re,c4Im);
  //printf("n=%d,re,im = %e,%e\n",nn,re,im);

  return CCTK_Cmplx(re,im);
}


void GetCoef_NBI_ScalarKerr_LeaverKerr(const CCTK_INT LeaverOrderRadial,
                                       const CCTK_COMPLEX lambda,
                                       const CCTK_INT mm,
                                       const CCTK_COMPLEX omega,
                                       const CCTK_REAL mass,
                                       const CCTK_REAL MBH,
                                       const CCTK_REAL aBH,
                                       CCTK_COMPLEX *Coef) {

  Coef[0] = CCTK_Cmplx(1.0,0.0);
  {
    const CCTK_COMPLEX alp = get_alpha_NBI_ScalarKerr_coef_LeaverKerr(0,
                                                                      omega,
                                                                      mass,
                                                                      MBH,
                                                                      aBH,
                                                                      lambda,
                                                                      mm);
    const CCTK_COMPLEX bet = get_beta_NBI_ScalarKerr_coef_LeaverKerr(0,
                                                                     omega,
                                                                     mass,
                                                                     MBH,
                                                                     aBH,
                                                                     lambda,
                                                                     mm);
    Coef[1] = CCTK_CmplxMul(-1.0,CCTK_CmplxMul(CCTK_CmplxDiv(bet,alp),Coef[0]));
  }

  for(int nn=1;nn<LeaverOrderRadial-1;++nn) {
    const CCTK_COMPLEX alp = get_alpha_NBI_ScalarKerr_coef_LeaverKerr(nn,
                                                                      omega,
                                                                      mass,
                                                                      MBH,
                                                                      aBH,
                                                                      lambda,
                                                                      mm);
    const CCTK_COMPLEX bet = get_beta_NBI_ScalarKerr_coef_LeaverKerr(nn,
                                                                     omega,
                                                                     mass,
                                                                     MBH,
                                                                     aBH,
                                                                     lambda,
                                                                     mm);
    const CCTK_COMPLEX gam = get_gamma_NBI_ScalarKerr_coef_LeaverKerr(nn,
                                                                      omega,
                                                                      mass,
                                                                      MBH,
                                                                      aBH,
                                                                      lambda,
                                                                      mm);

    //printf("nn = %d,(alp,bet,gam) = (%e,%e,%e,%e,%e,%e)\n",
    //    CCTK_CmplxReal(alp),CCTK_CmplxImag(alp),
    //    CCTK_CmplxReal(bet),CCTK_CmplxImag(bet),
    //    CCTK_CmplxReal(gam),CCTK_CmplxImag(gam));
    //printf("nn = %d,(bet/alp,gam/alp) = (%e,%e,%e,%e)\n",
    //    CCTK_CmplxReal(CCTK_CmplxDiv(bet,alp)),CCTK_CmplxImag(CCTK_CmplxDiv(bet,alp)),
    //    CCTK_CmplxReal(CCTK_CmplxDiv(gam,alp)),CCTK_CmplxImag(CCTK_CmplxDiv(gam,alp)));

      Coef[nn+1] = CCTK_CmplxMul(CCTK_Cmplx(-1.0,0.0),
                      CCTK_CmplxAdd(CCTK_CmplxMul(CCTK_CmplxDiv(bet,alp),Coef[nn]),
                                    CCTK_CmplxMul(CCTK_CmplxDiv(gam,alp),Coef[nn-1])));
  }
  return;
}

CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverKerr(const CCTK_REAL rw,
                                                          const CCTK_INT LeaverOrderRadial,
                                                          const CCTK_COMPLEX omega,
                                                          const CCTK_REAL mass,
                                                          const CCTK_REAL MBH,
                                                          const CCTK_REAL aBH,
                                                          const CCTK_COMPLEX lambda,
                                                          const CCTK_INT mm,
                                                          const CCTK_COMPLEX *Coef)
{
  const CCTK_REAL rp = MBH + sqrt(MBH*MBH - aBH*aBH);
  const CCTK_REAL rm = MBH - sqrt(MBH*MBH - aBH*aBH);
  const CCTK_REAL omegac = aBH*mm/(2.0*MBH*aBH);
  const CCTK_COMPLEX sigma = CCTK_CmplxDiv(CCTK_CmplxMul(2.0*rp,CCTK_CmplxSub(omega,omegac)),
                                           rp-rm);
  const CCTK_COMPLEX omega2 = CCTK_CmplxMul(omega,omega);
  CCTK_COMPLEX qq = CCTK_CmplxSqrt(CCTK_CmplxSub(mass*mass,omega2));
  if(CCTK_CmplxReal(qq) > 0.0) {
    qq = CCTK_CmplxMul(-1.0,qq);
  }
  const CCTK_COMPLEX chi = CCTK_CmplxDiv(CCTK_CmplxSub(mass*mass,CCTK_CmplxMul(2.0,omega2)),qq);

  const CCTK_COMPLEX fac1 = CCTK_CmplxExp(CCTK_CmplxMul(sigma,
                                                CCTK_CmplxMul(CCTK_Cmplx(0.0,-1.0),log(fabs(rw-rp)))));
  const CCTK_COMPLEX fac2
    = CCTK_CmplxExp(CCTK_CmplxMul(log(fabs(rw-rm)),
                          CCTK_CmplxAdd(CCTK_CmplxMul(sigma,CCTK_Cmplx(1.0,0.0)),
                                        CCTK_CmplxSub(chi,1.0))));
  const CCTK_COMPLEX fac3 = CCTK_CmplxExp(CCTK_CmplxMul(qq,rw));

  //printf("fac1,fac2,fac3 = (%e,%e),(%e,%e),(%e,%e)\n",
  //    CCTK_CmplxReal(fac1),CCTK_CmplxImag(fac1),
  //    CCTK_CmplxReal(fac2),CCTK_CmplxImag(fac2),
  //    CCTK_CmplxReal(fac3),CCTK_CmplxImag(fac3));

  CCTK_COMPLEX sum = CCTK_Cmplx(0.0,0.0);
  for(int nn=0;nn<LeaverOrderRadial;++nn) {
    const CCTK_REAL uu = pow(fabs((rw - rp)/(rw - rm)),nn);
    sum = CCTK_CmplxAdd(sum,CCTK_CmplxMul(Coef[nn],uu));
    //printf("nn=%d,uu=%e,Coef = (%e,%e)\n",nn,uu,CCTK_CmplxReal(Coef[nn]),CCTK_CmplxImag(Coef[nn]));
    //printf("nn=%d,sum = (%e,%e)\n",nn,CCTK_CmplxReal(sum),CCTK_CmplxImag(sum));
  }

   return CCTK_CmplxMul(fac1,CCTK_CmplxMul(fac2,CCTK_CmplxMul(fac3,sum)));
}

CCTK_COMPLEX get_max_NBI_ScalarKerr_LeaverKerr(const CCTK_INT Nrmax,
                                               const CCTK_REAL rmin,
                                               const CCTK_REAL rmax,
                                               const CCTK_INT LeaverOrderRadial,
                                               const CCTK_COMPLEX omega,
                                               const CCTK_REAL mass,
                                               const CCTK_REAL MBH,
                                               const CCTK_REAL aBH,
                                               const CCTK_COMPLEX lambda,
                                               const CCTK_INT mm,
                                               const CCTK_COMPLEX *Coef) {
   const CCTK_REAL dr = (rmax - rmin)/(Nrmax-1);
   CCTK_COMPLEX phi_max = CCTK_Cmplx(0.0, 0.0);
   for(int nn=1;nn<=Nrmax;++nn)
   {
      const CCTK_REAL rw = nn*dr + rmin;
      const CCTK_COMPLEX phil = get_radial_profile_NBI_ScalarKerr_LeaverKerr(rw,
                                                                             LeaverOrderRadial,
                                                                             omega,
                                                                             mass,
                                                                             MBH,
                                                                             aBH,
                                                                             lambda,
                                                                             mm,
                                                                             Coef);
      //printf("nn=%d,rw = %e,phil = (%e,%e)\n",nn,rw,CCTK_CmplxReal(phil),CCTK_CmplxImag(phil));
      if(fabs(CCTK_CmplxReal(phil)) > fabs(CCTK_CmplxReal(phi_max)))
      {
         phi_max = phil;
      }
   }

   return phi_max;
}


