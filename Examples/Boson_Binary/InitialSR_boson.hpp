/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "CoordinateTransformations.hpp"
#include "MatterCCZ4.hpp"
#include "TensorAlgebra.hpp"
#include <fstream>
#include "SphericalHarmonics.hpp"
#include <cmath>
#include <complex>

//! Class which sets the initial scalar field matter config
class InitialSR_boson
{
  public:
    InitialSR_boson(const double a_L, const double a_dx,
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
        //double bh_mass;   //!< Mass on the initial BH
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of grid, for working out coords if neeeded

    };

    //! Function to compute the value of all the initial vars on the grid
    //template <class data_t> void compute(Cell<data_t> current_cell) const
    void compute(Cell<double> current_cell) const
    {

        // may want to add to existing values so load the vars
        ADMConformalVars::VarsWithGauge<double> vars;
        // assign them zero to be sure they are zeroed
        // including at boundaries where required
        VarsTools::assign(vars, 0.);
      
        // where am I?
        // Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        const Coordinates<double> coords(current_cell, m_dx, m_center);

        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        // data_t r = coords.get_radius();
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

        double phi_ang_real = real(spin_Y_l_m);

        double phi_tot = phi * phi_ang_real;


        // set the field values, as constant in space
        // ScalarField<Potential>::Vars<data_t> scalar_vars;

        // scalar_vars.phi = phi_tot;
        // scalar_vars.Pi = 0.0;

        // Set the background for a Schwazschild BH in isotropic coords
        //data_t psi = 1.0 + 0.5 * m_params.bh_mass / r;
        //data_t chi = pow(psi, -4.0);

        // calculate the appropriate value of K to solve the constraints
        //data_t V_of_phi, dVdphi;
        //m_potential.compute_potential(V_of_phi, dVdphi, scalar_vars);
        //data_t K_squared = 24.0 * M_PI * m_G_Newton * V_of_phi;
        //data_t K = sqrt(K_squared);

        // store the vars
        current_cell.store_vars(phi_tot, c_phi);
        current_cell.store_vars(0.0, c_Pi);
        //current_cell.store_vars(scalar_vars);
        //current_cell.store_vars(K, c_K);
        //current_cell.store_vars(chi, c_chi);
        //current_cell.store_vars(1.0, c_h11);
        //current_cell.store_vars(1.0, c_h22);
        //current_cell.store_vars(1.0, c_h33);
        //current_cell.store_vars(1.0, c_lapse);
    }

  protected:
    const double m_dx;
    // const double m_G_Newton;
    // const Potential m_potential;
    const double m_L;
    const double m_spacing;
    const std::array<double, CH_SPACEDIM> m_center;
    const std::vector<double> m_phi_values;
    const params_t m_params; //!< The matter initial condition params

};

#endif /* INITIALSCALARDATA_HPP_ */
