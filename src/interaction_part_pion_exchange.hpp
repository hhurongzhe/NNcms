#pragma once
#ifndef INTERACTION_PART_PE_HPP
#define INTERACTION_PART_PE_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"

namespace interaction_part_pion_exchange
{
    constexpr double PI = 3.141592653589793;

    double get_isospin_factor(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz)
    {
        double factor = 1.0;
        if (tz == 0 && ((l_initial + s) % 2 != 0))
        {
            factor = -3.0;
        }
        return factor;
    }

    // LO one-pion exchange potential.
    std::vector<double> potential_one_pion_exchange(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const double &x, const NN::NN_configs &configs)
    {
        std::vector<double> f(6, 0.0); // [f1,f2,f3,f4,f5,f6] vector.
        double regulator_power = configs.n_reg_one_pion_exchange;

        double ga = configs.axial_current_coupling_constant;
        double gaga = ga * ga;
        double fpi = configs.pion_decay_constant;
        double ff = fpi * fpi;
        double frefactor = -gaga / 4.0 / ff;

        double pmag = p_initial;
        double ppmag = p_final;
        double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);

        double f1 = 0.0;
        double f2 = 0.0;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;
        double f6 = 0.0;

        double f6_ope_neutral = frefactor / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x +
                                             configs.mass_pion_neutral * configs.mass_pion_neutral);
        double f6_ope_charged = frefactor / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x +
                                             configs.mass_pion_charged * configs.mass_pion_charged);

        if (tz == 0)
        {
            // for np channel, taking into account CIB effect.
            double f6_ope_I0 = -f6_ope_neutral - 2.0 * f6_ope_charged;
            double f6_ope_I1 = -f6_ope_neutral + 2.0 * f6_ope_charged;
            if ((l_initial + s) % 2 == 0)
            {
                // total isospin I=1.
                f6 = f6 + f6_ope_I1;
            }
            else
            {
                // total isospin I=0.
                f6 = f6 + f6_ope_I0;
            }
        }
        if (tz != 0)
        {
            // for nn and pp channel, it's simple.
            f6 = f6 + f6_ope_neutral;
        }

        f[0] = f1;
        f[1] = f2;
        f[2] = f3;
        f[3] = f4;
        f[4] = f5;
        f[5] = f6;
        for (auto &component : f)
        {
            component *= regulator;
        }
        return f;
    }

    // NLO two-pion exchange potential.
    std::vector<double> potential_two_pion_exchange_nlo(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const double &x, const NN::NN_configs &configs)
    {
        std::vector<double> f(6, 0.0); // [f1,f2,f3,f4,f5,f6] vector.
        double regulator_power = configs.n_reg_two_pion_exchange_nlo;

        double ga = configs.axial_current_coupling_constant;
        double gaga = ga * ga;
        double gagagaga = gaga * gaga;
        double fpi = configs.pion_decay_constant;
        double ff = fpi * fpi;
        double ffff = ff * ff;
        double mpi2 = pow(configs.mass_pion_averaged, 2);
        double mpi4 = pow(configs.mass_pion_averaged, 4);

        double pmag = p_initial;
        double ppmag = p_final;
        double q2 = ppmag * ppmag + pmag * pmag - 2.0 * ppmag * pmag * x;
        double qmag = sqrt(q2);
        double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
        double isospin_factor = get_isospin_factor(l_final, l_initial, s, j, tz);

        double f1 = 0.0;
        double f2 = 0.0;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;
        double f6 = 0.0;

        double f1_fac = interaction_aPWD::loop_function_L(qmag, configs) / (384.0 * PI * PI * ffff);
        double f1_part1 = 4.0 * mpi2 * (1.0 + 4.0 * gaga - 5.0 * gagagaga);
        double f1_part2 = qmag * qmag * (1.0 + 10.0 * gaga - 23.0 * gagagaga);
        double f1_part3 = -48.0 * gagagaga * mpi4 / (4.0 * mpi2 + qmag * qmag);
        f1 = isospin_factor * f1_fac * (f1_part1 + f1_part2 + f1_part3);
        f6 = -3.0 * gagagaga / (64.0 * PI * PI * ffff) * interaction_aPWD::loop_function_L(qmag, configs);
        f2 = -qmag * qmag * f6;

        f[0] = f1;
        f[1] = f2;
        f[2] = f3;
        f[3] = f4;
        f[4] = f5;
        f[5] = f6;
        for (auto &component : f)
        {
            component *= regulator;
        }
        return f;
    }

    // N2LO two-pion exchange potential.
    std::vector<double> potential_two_pion_exchange_n2lo(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const double &x, const NN::NN_configs &configs)
    {
        std::vector<double> f(6, 0.0); // [f1,f2,f3,f4,f5,f6] vector.
        double regulator_power = configs.n_reg_two_pion_exchange_n2lo;

        double ga = configs.axial_current_coupling_constant;
        double gaga = ga * ga;
        double fpi = configs.pion_decay_constant;
        double ff = fpi * fpi;
        double ffff = ff * ff;
        double mpi2 = pow(configs.mass_pion_averaged, 2);

        double pmag = p_initial;
        double ppmag = p_final;
        double q2 = ppmag * ppmag + pmag * pmag - 2.0 * ppmag * pmag * x;
        double qmag = sqrt(q2);
        double regulator = interaction_aPWD::regulator_function(pmag, ppmag, regulator_power, configs);
        double isospin_factor = get_isospin_factor(l_final, l_initial, s, j, tz);

        double f1 = 0.0;
        double f2 = 0.0;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;
        double f6 = 0.0;

        double f1_part1 = 3.0 * gaga / (16.0 * PI * ffff);
        double f1_part2 = 2.0 * mpi2 * (configs.c3 - 2.0 * configs.c1) + configs.c3 * qmag * qmag;
        double f1_part3 = 2.0 * mpi2 + qmag * qmag;
        double f1_part4 = interaction_aPWD::loop_function_A(qmag, configs);
        f1 = f1_part1 * f1_part2 * f1_part3 * f1_part4;
        f6 = -isospin_factor * gaga / (32.0 * PI * ffff) * configs.c4 * (4.0 * mpi2 + qmag * qmag) *
             interaction_aPWD::loop_function_A(qmag, configs);
        f2 = -qmag * qmag * f6;

        f[0] = f1;
        f[1] = f2;
        f[2] = f3;
        f[3] = f4;
        f[4] = f5;
        f[5] = f6;
        for (auto &component : f)
        {
            component *= regulator;
        }
        return f;
    }

} // end namespace interaction_part_pion_exchange

#endif // INTERACTION_PART_PE_HPP
