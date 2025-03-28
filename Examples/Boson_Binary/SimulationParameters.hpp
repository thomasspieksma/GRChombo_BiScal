/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialSR_boson.hpp"
#include "Potential.hpp"
#include "KerrBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("spheroidicity_param", initial_params.spheroidicity_param, 0.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 0.2);
        // pp.load("kerr_mass", kerr_params.mass);
        // pp.load("kerr_spin", kerr_params.spin);
        // pp.load("kerr_center", kerr_params.center, center);
        // pp.load("kerr_spin_direction", kerr_params.spin_direction,
                // {0., 0., 1.});
        // pp.load("initial_data_prefix", initial_data_prefix);
        pp.load("inner_r", inner_r, 0.0);
        pp.load("outer_r", outer_r, 250.0);

        // Lineout params
        pp.load("lineout_num_points", lineout_num_points, 10);

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * kerr_params.mass);

#endif
    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
    }

    // Initial data for matter and potential and BH
    double G_Newton, inner_r, outer_r;
    int lineout_num_points;
    InitialSR_boson::params_t initial_params;
    KerrBH::params_t kerr_params;
    Potential::params_t potential_params;
    std::string initial_data_prefix;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
