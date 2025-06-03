#pragma once
#ifndef INTERACTION_PART_CONTACT_HPP
#define INTERACTION_PART_CONTACT_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"

namespace interaction_part_contact
{

    // LO contact potential.
    double potential_contact_lo(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const NN::NN_configs &configs)
    {
        double pmag = p_initial;
        double ppmag = p_final;

        if (l_final == 0 && l_initial == 0 && s == 0 && j == 0) // 1S0 channel, there is CIB.
        {
            double regulator_power = configs.n_reg_Ctilde_1s0;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);

            if (tz == -1)
            {
                return regulator * configs.Ctilde_1s0_pp;
            }
            if (tz == 0)
            {
                return regulator * configs.Ctilde_1s0_np;
            }
            if (tz == 1)
            {
                return regulator * configs.Ctilde_1s0_nn;
            }
        }
        else if (l_final == 0 && l_initial == 0 && s == 1 && j == 1) // 3S1 channel.
        {
            double regulator_power = configs.n_reg_Ctilde_3s1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);

            return regulator * configs.Ctilde_3s1;
        }
        else
        {
            return 0.0;
        }
        return 0.0;
    }

    // NLO contact potential.
    double potential_contact_nlo(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const NN::NN_configs &configs)
    {
        double pmag = p_initial;
        double ppmag = p_final;

        if (l_final == 0 && l_initial == 0 && s == 0 && j == 0) // 1S0 channel.
        {
            double regulator_power = configs.n_reg_C_1s0;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_1s0 * (pmag * pmag + ppmag * ppmag);
        }
        else if (l_final == 1 && l_initial == 1 && s == 1 && j == 0) // 3P0 channel.
        {
            double regulator_power = configs.n_reg_C_3p0;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3p0 * pmag * ppmag;
        }
        else if (l_final == 1 && l_initial == 1 && s == 0 && j == 1) // 1P1 channel.
        {
            double regulator_power = configs.n_reg_C_1p1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_1p1 * pmag * ppmag;
        }
        else if (l_final == 1 && l_initial == 1 && s == 1 && j == 1) // 3P1 channel.
        {
            double regulator_power = configs.n_reg_C_3p1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3p1 * pmag * ppmag;
        }
        else if (l_final == 0 && l_initial == 0 && s == 1 && j == 1) // 3S1 channel.
        {
            double regulator_power = configs.n_reg_C_3s1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3s1 * (pmag * pmag + ppmag * ppmag);
        }
        else if (l_final == 0 && l_initial == 2 && s == 1 && j == 1) // 3S1-3D1 channel.
        {
            double regulator_power = configs.n_reg_C_3sd1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3sd1 * pmag * pmag;
        }
        else if (l_final == 2 && l_initial == 0 && s == 1 && j == 1) // 3D1-3S1 channel.
        {
            double regulator_power = configs.n_reg_C_3sd1;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3sd1 * ppmag * ppmag;
        }
        else if (l_final == 1 && l_initial == 1 && s == 1 && j == 2) // 3P2 channel.
        {
            double regulator_power = configs.n_reg_C_3p2;
            double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
            return regulator * configs.C_3p2 * pmag * ppmag;
        }
        else
        {
            return 0.0;
        }
    }

} // end namespace interaction_part_contact

#endif // INTERACTION_PART_CONTACT_HPP
