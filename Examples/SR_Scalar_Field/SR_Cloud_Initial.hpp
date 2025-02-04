/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SRCLOUDINITIAL_HPP_
#define SRCLOUDINITIAL_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <fstream>
#include "SphericalHarmonics.hpp"
#include <cmath>
#include <complex>
// Class which creates initial conditions for a scalar boson cloud formed through superradiance
// We make use of an external data-file with radial data
// See Appendix A of arXiv:2306.16447

class SR_Cloud_Initial
{
  public:
    SR_Cloud_Initial(const double a_L, const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const double spacing,
                      const std::vector<double> phi_values)
        : m_L(a_L), m_dx(a_dx), m_spacing(spacing), m_center(a_center),
          m_phi_values(phi_values)
    
    {
    }
    struct params_t
        {
        double spheroidicity_param; 
    };

    //! Function to compute the value of all the initial vars on the grid
    void compute(Cell<double> current_cell) const
    {
        // may want to add to existing values so load the vars
        // ADMConformalVars::VarsWithGauge<double> vars;
        // assign them zero to be sure they are zeroed
        // including at boundaries where required
        // VarsTools::assign(vars, 0.);

        // Define Coordinates
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r = coords.get_radius();
        double rho2 = max(x * x + y * y, 1e-12);
        double r2sin2theta = rho2;

        // Interpolate data from read in values
        const int indxL = static_cast<int>(floor(r / m_spacing));
        const int indxH = static_cast<int>(ceil(r / m_spacing));

        // the field values
        double phi =
        m_phi_values[indxL] +
            (r / m_spacing - indxL) * (m_phi_values[indxH] - m_phi_values[indxL]);
        
        int ell_SH = 1;
        int m_SH = 1;
        int s_SH = 0;
        
        // SphericalHarmonics::Y_lm_t<double> spi_Y_Min2Plusl_m = SphericalHarmonics::spin_Y_lm(x, y, z, s_SH, -2.0 + ell_SH, m_SH);
        SphericalHarmonics::Y_lm_t<double> spi_Y_2Plusl_m = SphericalHarmonics::spin_Y_lm(x, y, z, s_SH, 2.0 + ell_SH, m_SH);
        // SphericalHarmonics::Y_lm_t<double> spi_Y_Min4Plusl_m = SphericalHarmonics::spin_Y_lm(x, y, z, s_SH, -4.0 + ell_SH, m_SH);
        SphericalHarmonics::Y_lm_t<double> spi_Y_l_m = SphericalHarmonics::spin_Y_lm(x, y, z, s_SH, ell_SH, m_SH);
        SphericalHarmonics::Y_lm_t<double> spi_Y_4Plusl_m = SphericalHarmonics::spin_Y_lm(x, y, z, s_SH, 4.0 + ell_SH, m_SH);

        const   complex<double> i(0.0,1.0);
        complex<double> spin_Y_Min2Plusl_m = 0.0; //spi_Y_Min2Plusl_m.Real + i * spi_Y_Min2Plusl_m.Im;
        complex<double> spin_Y_2Plusl_m = spi_Y_2Plusl_m.Real + i * spi_Y_2Plusl_m.Im;
        complex<double> spin_Y_Min4Plusl_m = 0.0; //spi_Y_Min4Plusl_m.Real + i * spi_Y_Min4Plusl_m.Im;
        complex<double> spin_Y_l_m = spi_Y_l_m.Real + i * spi_Y_l_m.Im;
        complex<double> spin_Y_4Plusl_m = spi_Y_4Plusl_m.Real + i * spi_Y_4Plusl_m.Im;


        complex<double> nom_gamm_sq = -pow(3.0 + 2.0 * ell_SH, 2.0) * sqrt(5.0 + 2.0 * ell_SH) * sqrt(ell_SH * ell_SH - m_SH * m_SH) * sqrt(1.0 - 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * spin_Y_Min2Plusl_m + pow(1.0 - 2.0 * ell_SH, 2.0) * sqrt(-3.0 + 2.0 * ell_SH) * sqrt(1.0 + 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(4.0 + 4.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * spin_Y_2Plusl_m;
        complex<double> denom_gamm_sq = 2.0 * sqrt(-3.0 + 2.0 * ell_SH) * sqrt(1.0 + 2.0 * ell_SH) * sqrt(5.0 + 2.0 * ell_SH) * pow(-3.0 + 4.0 * ell_SH + 4.0 * pow(ell_SH, 2.0), 2.0);
        complex<double> gamm_qaur_term1_nom = sqrt(ell_SH * ell_SH - m_SH * m_SH) * sqrt(9.0 - 6.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(4.0 - 4.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(1.0 - 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * spin_Y_Min4Plusl_m;
        complex<double> gamm_qaur_term1_denom = sqrt(-7.0 + 2.0 * ell_SH) * (-5.0 + 2.0 * ell_SH) * sqrt(1.0 + 2.0 * ell_SH) * pow(3.0 -8.0 * ell_SH + 4.0 * ell_SH * ell_SH, 2.0);
        complex<double> gamm_qaur_term1 = gamm_qaur_term1_nom / gamm_qaur_term1_denom;
        
        complex<double> gamm_qaur_term2_nom = 8.0 * sqrt(ell_SH * ell_SH - m_SH * m_SH) * sqrt(1.0 - 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * (-1.0 + 4.0 * m_SH * m_SH) * spin_Y_Min2Plusl_m;
        complex<double> gamm_qaur_term2_denom = pow(1.0 - 2.0 * ell_SH, 4.0) * (-5.0 + 2.0 * ell_SH) * sqrt(-3.0 + 2.0 * ell_SH) * sqrt(1.0 + 2.0 * ell_SH) * (3.0 + 2.0 * ell_SH);
        complex<double> gamm_qaur_term2 = gamm_qaur_term2_nom / gamm_qaur_term2_denom;
        
        complex<double> gamm_qaur_term3_nom = ((ell_SH * ell_SH - m_SH * m_SH) * (1.0 - 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) / (pow(1.0 - 2.0 * ell_SH, 4.0) * (-3.0 + 2.0 * ell_SH)) + ((1.0 + 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * (4.0 + 4.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH)) / (pow(3.0 + 2.0 * ell_SH, 4.0) * (5.0 + 2.0 * ell_SH))) * spin_Y_l_m;
        complex<double> gamm_qaur_term3_denom = 1.0 + 2.0 * ell_SH;
        complex<double> gamm_qaur_term3 = gamm_qaur_term3_nom / gamm_qaur_term3_denom;
        
        complex<double> gamm_qaur_term4_nom = 8.0 * sqrt(1.0 + 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(4.0 + 4.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * (-1.0 + 4.0 * m_SH * m_SH) * spin_Y_2Plusl_m;
        complex<double> gamm_qaur_term4_denom = (-1.0 + 2.0 * ell_SH) * sqrt(1.0 + 2.0 * ell_SH) * pow(3.0 + 2.0 * ell_SH, 4.0) * sqrt(5.0 + 2.0 * ell_SH) * (7.0 + 2.0 * ell_SH);
        complex<double> gamm_qaur_term4 = gamm_qaur_term4_nom / gamm_qaur_term4_denom;
        
        complex<double> gamm_qaur_term5_nom = sqrt(1.0 + 2.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(4.0 + 4.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(9.0 + 6.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * sqrt(16.0 + 8.0 * ell_SH + ell_SH * ell_SH - m_SH * m_SH) * spin_Y_4Plusl_m;
        complex<double> gamm_qaur_term5_denom = sqrt(1.0 + 2.0 * ell_SH) * pow(3.0 + 2.0 * ell_SH, 2.0) * pow(5.0 + 2.0 * ell_SH, 2.0) * (7.0 + 2.0 * ell_SH);
        complex<double> gamm_qaur_term5 = gamm_qaur_term5_nom / gamm_qaur_term5_denom;
        
        complex<double> phi_ang = spin_Y_l_m + m_params.spheroidicity_param * nom_gamm_sq / denom_gamm_sq + 1.0/8.0 * m_params.spheroidicity_param * m_params.spheroidicity_param * (gamm_qaur_term1 - gamm_qaur_term2 - gamm_qaur_term3 + gamm_qaur_term4 + gamm_qaur_term5);

        double phi_ang_real = real(phi_ang);

        double phi_tot = phi * phi_ang_real;

        // store the metric vars
        // current_cell.store_vars(vars);
        // store the matter vars
        current_cell.store_vars(phi_tot, c_phi);
        current_cell.store_vars(0.0, c_Pi);
    }

  protected:
    const double m_dx;
    const double m_L;
    const double m_spacing;
    const params_t m_params;
    const std::array<double, CH_SPACEDIM> m_center;
    const std::vector<double> m_phi_values;
};

#endif /* SRCLOUDINITIAL_HPP_ */
