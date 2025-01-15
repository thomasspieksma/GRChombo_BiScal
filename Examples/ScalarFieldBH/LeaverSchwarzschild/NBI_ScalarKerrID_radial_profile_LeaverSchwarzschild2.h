#include <math.h>
#include <stdbool.h>

#include "NBI_ScalarKerrID_coef_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_LeaverSolver_LeaverSchwarzschild.h"
#include "NBI_ScalarKerrID_radial_profile_LeaverSchwarzschild.h"

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
                                                     CCTK_INT N_grid_search_max);
CCTK_COMPLEX get_radial_profile_NBI_ScalarKerr_LeaverSchwarzschild2(CCTK_REAL rw,
                                                                    CCTK_INT NmaxLeaver,
                                                                    CCTK_COMPLEX omega,
                                                                    CCTK_REAL mass,
                                                                    CCTK_REAL MBH,
                                                                    CCTK_REAL ll,
                                                                    CCTK_COMPLEX *Coef);

