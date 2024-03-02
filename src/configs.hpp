#pragma once
#ifndef CONFIGS_HPP
#define CONFIGS_HPP

#include "inifile.hpp"
#include "basic_math.hpp"
#include "constants.hpp"

namespace NN
{

    struct NN_configs
    {

        // ***** interaction section *****
        double axial_current_coupling_constant;
        double pion_decay_constant;
        double c1, c3, c4;
        double Ctilde_1s0_pp, Ctilde_1s0_nn, Ctilde_1s0_np, Ctilde_3s1;
        double C_1s0, C_3s1, C_1p1, C_3p0, C_3p1, C_3sd1, C_3p2;
        double Lambda, Lambda_tilde;

        int n_reg_Ctilde_1s0, n_reg_Ctilde_3s1, n_reg_C_1s0, n_reg_C_3s1, n_reg_C_1p1, n_reg_C_3p0, n_reg_C_3p1, n_reg_C_3sd1, n_reg_C_3p2;
        int n_reg_one_pion_exchange, n_reg_two_pion_exchange_nlo, n_reg_two_pion_exchange_n2lo;

        // ***** meson masses section ****
        double mass_pion_charged;
        double mass_pion_neutral;
        double mass_pion_averaged;

        // ***** baryon masses section ****
        double mass_proton;
        double mass_neutron;
        double mass_nucleon;

        // ***** numerical parameters section ****

        // number of momentum mesh points.
        int mesh_points_number;
        std::vector<double> momentum_mesh_points;
        std::vector<double> momentum_mesh_weights;

        // angular mesh points, we fix them.
        int angular_mesh_number = 24;
        std::vector<double> angular_mesh_points;
        std::vector<double> angular_mesh_weights;

        double k_max;
        int J_max;
        bool use_omp;
        int n_omp_threads;

        // ***** output section *****
        std::string result_dir;
        std::string result_name;

        // constructor
        NN_configs(const inifile_system::inifile &ini);

        // generate a result file name
        std::string result_file() const;
    };

    NN_configs::NN_configs(const inifile_system::inifile &ini)
    {

        // ***** interaction section *****
        auto sec = ini.section("interaction");
        axial_current_coupling_constant = sec.get_double("axial_current_coupling_constant");
        pion_decay_constant = sec.get_double("pion_decay_constant");

        c1 = sec.get_double("c1") * 1e-3;
        c3 = sec.get_double("c3") * 1e-3;
        c4 = sec.get_double("c4") * 1e-3;

        Ctilde_1s0_pp = sec.get_double("Ctilde_1s0_pp") * 1e-2;
        Ctilde_1s0_nn = sec.get_double("Ctilde_1s0_nn") * 1e-2;
        Ctilde_1s0_np = sec.get_double("Ctilde_1s0_np") * 1e-2;
        Ctilde_3s1 = sec.get_double("Ctilde_3s1") * 1e-2;

        C_1s0 = sec.get_double("C_1s0") * 1e-8;
        C_3s1 = sec.get_double("C_3s1") * 1e-8;
        C_1p1 = sec.get_double("C_1p1") * 1e-8;
        C_3p0 = sec.get_double("C_3p0") * 1e-8;
        C_3p1 = sec.get_double("C_3p1") * 1e-8;
        C_3sd1 = sec.get_double("C_3sd1") * 1e-8;
        C_3p2 = sec.get_double("C_3p2") * 1e-8;

        Lambda = sec.get_double("Lambda");
        Lambda_tilde = sec.get_double("Lambda_tilde");

        n_reg_Ctilde_1s0 = sec.get_int("n_reg_Ctilde_1s0");
        n_reg_Ctilde_3s1 = sec.get_int("n_reg_Ctilde_3s1");
        n_reg_C_1s0 = sec.get_int("n_reg_C_1s0");
        n_reg_C_3s1 = sec.get_int("n_reg_C_3s1");
        n_reg_C_1p1 = sec.get_int("n_reg_C_1p1");
        n_reg_C_3p0 = sec.get_int("n_reg_C_3p0");
        n_reg_C_3p1 = sec.get_int("n_reg_C_3p1");
        n_reg_C_3sd1 = sec.get_int("n_reg_C_3sd1");
        n_reg_C_3p2 = sec.get_int("n_reg_C_3p2");

        n_reg_one_pion_exchange = sec.get_int("n_reg_one_pion_exchange");
        n_reg_two_pion_exchange_nlo = sec.get_int("n_reg_two_pion_exchange_nlo");
        n_reg_two_pion_exchange_n2lo = sec.get_int("n_reg_two_pion_exchange_n2lo");

        // ***** meson masses section ****
        sec = ini.section("meson-masses");
        mass_pion_charged = sec.get_double("mass_pion_charged");
        mass_pion_neutral = sec.get_double("mass_pion_neutral");
        mass_pion_averaged = sec.get_double("mass_pion_averaged");

        // ***** baryon masses section ****
        sec = ini.section("baryon-masses");
        mass_proton = sec.get_double("mass_proton");
        mass_neutron = sec.get_double("mass_neutron");
        mass_nucleon = sec.get_double("mass_nucleon");

        // ***** numerical parameters section ****
        sec = ini.section("numerical-parameters");
        mesh_points_number = sec.get_int("mesh_points_number");
        k_max = constants::hbarc * sec.get_double("k_max");
        J_max = sec.get_int("J_max");
        use_omp = sec.get_bool("use_omp");
        n_omp_threads = sec.get_int("n_omp_threads");
        if (!use_omp)
        {
            n_omp_threads = 1;
        }
        momentum_mesh_points = basic_math::gauss_legendre_nodes_interval(0, k_max, mesh_points_number);
        momentum_mesh_weights = basic_math::gauss_legendre_weights_interval(0, k_max, mesh_points_number);
        angular_mesh_points = basic_math::gauss_legendre_nodes(angular_mesh_number);
        angular_mesh_weights = basic_math::gauss_legendre_weights(angular_mesh_number);

        // ***** output section *****
        sec = ini.section("output");
        result_dir = sec.get_string("result_dir");
        result_name = sec.get_string("result_name");
        if (result_dir.back() != '/')
        {
            result_dir += "/";
        }
    };

    std::string file_stem(const std::string &file)
    {
        auto p1 = file.find_last_of('/') + 1;
        auto p2 = file.find_last_of('.');
        return file.substr(p1, p2 - p1);
    }

    std::string NN_configs::result_file() const
    {
        std::ostringstream oss;
        oss << result_dir << result_name << ".dat";
        return oss.str();
    }

} // namespace NN

#endif // CONFIGS_HPP