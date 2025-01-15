///////////////////////////////////////////////
#ifndef INCLUDE_HEAD_NBI_SCALARKERRID_SPHEROIDAL_HARMONICS2
#define INCLUDE_HEAD_NBI_SCALARKERRID_SPHEROIDAL_HARMONICS2
//////////////////////////////////////////////


#include <math.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_alpha(const CCTK_INT nn,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm,
                                                      const CCTK_COMPLEX gg);

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_beta(const CCTK_INT nn,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm,
                                                      const CCTK_COMPLEX gg);

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_gamma(const CCTK_INT nn,
                                                      const CCTK_COMPLEX lambda,
                                                      const CCTK_INT mm,
                                                      const CCTK_COMPLEX gg);

void nbi_scalarkerr_SpheroidalHarmonics_GetCoef(const CCTK_INT LeaverOrderAng,
                                                const CCTK_INT lambda,
                                                const CCTK_INT mm,
                                                const CCTK_COMPLEX gg,
                                                CCTK_COMPLEX *Coef);

CCTK_COMPLEX nbi_scalarkerr_SpheroidalHarmonics_GetValue(const CCTK_INT LeaverOrderAng,
                                                         const CCTK_INT lambda,
                                                         const CCTK_INT mm,
                                                         const CCTK_COMPLEX gg,
                                                         CCTK_COMPLEX *Coef,
                                                         const CCTK_REAL th,
                                                         const CCTK_REAL ph);
#endif
