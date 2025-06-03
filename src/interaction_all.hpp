#pragma once
#ifndef ALL_INTERACTION_HPP
#define ALL_INTERACTION_HPP

#include "interaction_part_pion_exchange.hpp"
#include "interaction_part_contact.hpp"
#include "lib_define.hpp"
#include <omp.h>

namespace interaction_all
{
    constexpr double twopicubic = 248.0502134423985614038105; // (2*Pi)^3

    double potential_chiral(const int &l_final, const int &l_initial, const int &s, const int &j, const int &tz, const double &p_final, const double &p_initial, const NN::NN_configs &configs)
    {
        double nucleon_mass = configs.mass_nucleon;
        double e_final = sqrt(nucleon_mass * nucleon_mass + p_final * p_final);
        double e_initial = sqrt(nucleon_mass * nucleon_mass + p_initial * p_initial);
        double relativity_factor = nucleon_mass / sqrt(e_final * e_initial);
        double temp = 0.0;

        // contact terms, already partial-wave projected.
        // lo terms.
        temp = temp + interaction_part_contact::potential_contact_lo(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
        // nlo terms.
        temp = temp + interaction_part_contact::potential_contact_nlo(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
        // n2lo terms.
        // there is no n2lo contact terms.

        std::vector<double> f_component_vec(6, 0.0); // [f1,f2,f3,f4,f5,f6] vector.
        double x, w, fa;
        std::vector<double> one_pion_exchange;
        std::vector<double> two_pion_exchange_nlo;
        std::vector<double> two_pion_exchange_n2lo;

#pragma omp parallel for private(f_component_vec, x, w, fa, one_pion_exchange, two_pion_exchange_nlo, two_pion_exchange_n2lo) reduction(+ : temp) schedule(dynamic)
        // pion-exchange terms, need to do PWD.
        for (size_t idx_angle = 0; idx_angle < configs.angular_mesh_number; idx_angle = idx_angle + 1)
        {
            f_component_vec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            x = configs.angular_mesh_points[idx_angle]; // x=cos(theta), where theta is the angle between p_final and p_initial.
            w = configs.angular_mesh_weights[idx_angle];

            // LO one-pion-exchange term.
            one_pion_exchange = interaction_part_pion_exchange::potential_one_pion_exchange(l_final, l_initial, s, j, tz, p_final, p_initial, x, configs);
            // NLO two-pion-exchange term.
            two_pion_exchange_nlo = interaction_part_pion_exchange::potential_two_pion_exchange_nlo(l_final, l_initial, s, j, tz, p_final, p_initial, x, configs);
            // N2LO two-pion-exchange term.
            two_pion_exchange_n2lo = interaction_part_pion_exchange::potential_two_pion_exchange_n2lo(l_final, l_initial, s, j, tz, p_final, p_initial, x, configs);

            for (size_t idx_f = 0; idx_f < f_component_vec.size(); idx_f = idx_f + 1)
            {
                // adding lo terms.
                f_component_vec[idx_f] += one_pion_exchange[idx_f];
                // adding nlo terms.
                f_component_vec[idx_f] += two_pion_exchange_nlo[idx_f];
                // adding n2lo terms.
                f_component_vec[idx_f] += two_pion_exchange_n2lo[idx_f];
            }
            // perform aPWD after adding up all terms. (independent of tz)
            fa = interaction_aPWD::potential_auto(l_final, l_initial, s, j, p_final, p_initial, x, f_component_vec);
            temp = temp + fa * w;
        }

        // apply a relativity-factor and a normalization constant (2Pi)^3.
        temp = temp * relativity_factor / twopicubic;
        return temp;
    }

} // namespace interaction_all

#endif // ALL_INTERACTION_HPP
