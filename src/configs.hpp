#pragma once
#ifndef CONFIGS_HPP
#define CONFIGS_HPP

#include "inifile.hpp"
#include "basic_math.hpp"
#include "constants.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

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
        int angular_mesh_number;
        std::vector<double> angular_mesh_points;
        std::vector<double> angular_mesh_weights;

        bool is_kernel;
        bool is_aside;
        bool use_omp;
        int n_omp_threads;

        std::vector<double> tlabs;

        std::vector<std::vector<int>> partial_waves;

        // ***** output section *****
        std::string result_dir;
        std::string result_name;

        // constructor
        NN_configs(const inifile_system::inifile &ini);

        // generate a result file name
        std::string result_file() const;

        // read momentum mesh from file
        void read_momentum_mesh(std::string file_momentum_mesh);

        // read tlabs table from file
        void read_table_tlabs(std::string file_tlabs);

        // read partial-waves from file
        void read_partial_waves(std::string file_pw);

        // get c.m. momentum coressponding to tlab
        double get_rel_mom(double tlab, int tz) const;
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
        is_kernel = sec.get_bool("is_kernel");
        is_aside = sec.get_bool("is_aside");
        use_omp = sec.get_bool("use_omp");
        n_omp_threads = sec.get_int("n_omp_threads");
        angular_mesh_number = sec.get_int("angular_mesh_number");
        if (!use_omp)
        {
            n_omp_threads = 1;
        }
        // set up angular mesh.
        angular_mesh_points = basic_math::gauss_legendre_nodes(angular_mesh_number);
        angular_mesh_weights = basic_math::gauss_legendre_weights(angular_mesh_number);

        // read momentum mesh.
        std::string file_momentum_mesh = "table_momentum_mesh.txt";
        read_momentum_mesh(file_momentum_mesh);

        // read tlabs if aside.
        std::string file_tlabs = "table_tlabs.txt";
        if (is_aside)
        {
            read_table_tlabs(file_tlabs);
        }

        // read partial-waves.
        std::string file_pw = "table_partial_waves.txt";
        read_partial_waves(file_pw);

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

    void NN_configs::read_momentum_mesh(std::string file_momentum_mesh)
    {
        std::ifstream file(file_momentum_mesh);
        if (!file.is_open())
        {
            std::cerr << "Error file_momentum_mesh: " << file_momentum_mesh << std::endl;
        }
        std::string line;
        // Skip comments
        while (std::getline(file, line) && line[0] == '#')
        {
        }
        // Read mesh_points_number
        std::istringstream iss(line);
        iss >> mesh_points_number;
        // Read momentum mesh points and weights
        double point, weight;
        while (file >> point >> weight)
        {
            momentum_mesh_points.push_back(point);
            momentum_mesh_weights.push_back(weight);
        }
        // Check if the number of points read matches mesh_points_number
        if (momentum_mesh_points.size() != static_cast<size_t>(mesh_points_number) ||
            momentum_mesh_weights.size() != static_cast<size_t>(mesh_points_number))
        {
            std::cerr << "Error: Number of points read doesn't match mesh_points_number." << std::endl;
        }
        file.close();
    }

    void NN_configs::read_table_tlabs(std::string file_tlabs)
    {
        std::ifstream file(file_tlabs);
        if (!file.is_open())
        {
            std::cerr << "Error file_tlabs: " << file_tlabs << std::endl;
        }
        std::string line;
        // Skip comments
        while (std::getline(file, line) && line[0] == '#')
        {
        }
        // Read tlabs
        double tlab;
        std::istringstream iss(line);
        iss >> tlab;
        tlabs.push_back(tlab);
        while (file >> tlab)
        {
            tlabs.push_back(tlab);
        }
        file.close();
    }

    void NN_configs::read_partial_waves(std::string file_pw)
    {
        std::ifstream file(file_pw);
        if (!file.is_open())
        {
            std::cerr << "Error file_pw: " << file_pw << std::endl;
        }
        std::string line;
        // Skip comments
        while (std::getline(file, line) && line[0] == '#')
        {
        }
        // Read partial-waves
        int lp, l, s, j, tz;
        std::vector<int> temp;
        std::istringstream iss(line);
        iss >> lp >> l >> s >> j >> tz;
        temp = {lp, l, s, j, tz};
        partial_waves.push_back(temp);
        while (file >> lp >> l >> s >> j >> tz)
        {
            temp = {lp, l, s, j, tz};
            partial_waves.push_back(temp);
        }
        file.close();
    }

    double NN_configs::get_rel_mom(double tlab, int tz) const
    {
        double q2; // square of c.m. momentum
        if (tz == -1)
        {
            q2 = 0.5 * mass_proton * tlab;
        }
        else if (tz == 0)
        {
            q2 = mass_proton * mass_proton * tlab * (tlab + 2.0 * mass_neutron) / ((mass_proton + mass_neutron) * (mass_proton + mass_neutron) + 2.0 * tlab * mass_proton);
        }
        else if (tz == 1)
        {
            q2 = 0.5 * mass_neutron * tlab;
        }
        else
        {
            std::cerr << "unknown isospin projection!" << std::endl;
        }
        double q;
        if (q2 < 0)
        {
            std::cout << "q2 < 0 !" << std::endl;
            q = std::sqrt(-q2);
        }
        else
        {
            q = std::sqrt(q2);
        }
        return q;
    }

} // namespace NN

#endif // CONFIGS_HPP