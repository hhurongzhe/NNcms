#include "lib_define.hpp"
#include "interaction_all.hpp"

void write_dat_single_channel(std::vector<int> this_channel, const NN::NN_configs &configs)
{
    int l_final, l_initial, s, j, tz;
    l_final = this_channel[0];
    l_initial = this_channel[1];
    s = this_channel[2];
    j = this_channel[3];
    tz = this_channel[4];
    std::string tz_name;
    if (tz == -1)
    {
        tz_name = "pp";
    }
    else if (tz == 0)
    {
        tz_name = "np";
    }
    else if (tz == 1)
    {
        tz_name = "nn";
    }
    else
    {
        tz_name = "???";
    }

    if (configs.is_kernel)
    {
        std::ostringstream oss;
        oss << configs.result_dir << configs.result_name << "-" << l_final << "-" << l_initial << "-" << s << "-" << j << "-" << tz_name << "-kernel.dat";
        auto file_name_this_channel = oss.str(); // file name for this partial-wave channel.
        std::cout << "writing: " << file_name_this_channel << std::endl;
        std::ofstream fp(file_name_this_channel);
        for (int idx_mom_bra = 0; idx_mom_bra < configs.mesh_points_number; idx_mom_bra = idx_mom_bra + 1)
        {
            for (int idx_mom_ket = 0; idx_mom_ket < configs.mesh_points_number; idx_mom_ket = idx_mom_ket + 1)
            {
                double p_final = configs.momentum_mesh_points[idx_mom_bra];
                double p_initial = configs.momentum_mesh_points[idx_mom_ket];
                double v_value = interaction_all::potential_chiral(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
                fp << std::fixed << " " << std::scientific << std::setprecision(17) << v_value;
            }
            fp << "\n";
        }
        fp.close();
    }
    if (configs.is_aside)
    {
        std::ostringstream oss;
        oss << configs.result_dir << configs.result_name << "-" << l_final << "-" << l_initial << "-" << s << "-" << j << "-" << tz_name << "-aside.dat";
        auto file_name_this_channel = oss.str(); // file name for this partial-wave channel.
        std::cout << "writing: " << file_name_this_channel << std::endl;
        std::ofstream fp(file_name_this_channel);
        for (int idx_tlab = 0; idx_tlab < configs.tlabs.size(); idx_tlab = idx_tlab + 1)
        {
            double tlab = configs.tlabs[idx_tlab];
            double p_final;
            double p_initial = configs.get_rel_mom(tlab, tz);
            double v_value;
            for (int idx_mom_bra = 0; idx_mom_bra < configs.mesh_points_number; idx_mom_bra = idx_mom_bra + 1)
            {
                p_final = configs.momentum_mesh_points[idx_mom_bra];
                v_value = interaction_all::potential_chiral(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
                fp << std::fixed << " " << std::scientific << std::setprecision(17) << v_value;
            }
            p_final = p_initial;
            v_value = interaction_all::potential_chiral(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
            fp << std::fixed << " " << std::scientific << std::setprecision(17) << v_value;
            fp << "\n";
        }
        fp.close();
    }
}

// write results in the output file.
void write_dat(const NN::NN_configs &configs)
{
    // write momentum mesh
    std::ostringstream oss_mom_mesh;
    oss_mom_mesh << configs.result_dir << configs.result_name << "-momentum-mesh.dat";
    auto file_mom_mesh = oss_mom_mesh.str();
    std::ofstream fp_mom_mesh(file_mom_mesh);
    fp_mom_mesh << "# momentum mesh points number;\n";
    fp_mom_mesh << "# momentum mesh points and weights;\n";
    fp_mom_mesh << configs.mesh_points_number << "\n";
    for (int i = 0; i < configs.momentum_mesh_points.size(); i = i + 1)
    {
        fp_mom_mesh << std::fixed << std::scientific << std::setprecision(17) << configs.momentum_mesh_points[i] << "\t" << std::setw(17) << configs.momentum_mesh_weights[i] << std::endl;
    }
    fp_mom_mesh.close();

    // write tlabs
    std::ostringstream oss_tlabs;
    oss_tlabs << configs.result_dir << configs.result_name << "-tlabs.dat";
    auto file_tlabs = oss_tlabs.str();
    std::ofstream fp_tlabs(file_tlabs);
    fp_tlabs << "# tlabs;\n";
    for (int i = 0; i < configs.tlabs.size(); i = i + 1)
    {
        fp_tlabs << std::fixed << std::scientific << std::setprecision(17) << configs.tlabs[i] << std::endl;
    }
    fp_tlabs.close();

    // TODO: write partial-waves
    std::ostringstream oss_pws;
    oss_pws << configs.result_dir << configs.result_name << "-partial-waves.dat";
    auto file_pws = oss_pws.str();
    std::ofstream fp_pws(file_pws);
    fp_pws << "# partial-waves: l' l s j tz;\n";
    for (int i = 0; i < configs.partial_waves.size(); i = i + 1)
    {
        fp_pws << configs.partial_waves[i][0] << " " << configs.partial_waves[i][1] << " " << configs.partial_waves[i][2] << " " << configs.partial_waves[i][3] << " " << configs.partial_waves[i][4] << "\n";
    }
    fp_pws.close();

    // generate channels.
    auto channels = configs.partial_waves;
    for (int idx_channel = 0; idx_channel < channels.size(); idx_channel = idx_channel + 1)
    {
        // write matrix elements for each channel.
        auto this_channel = channels[idx_channel];
        write_dat_single_channel(this_channel, configs);
    }
}

int main()
{
    // config file initializing.
    auto ini = inifile_system::inifile("inifile-cms.ini");
    if (!ini.good())
    {
        std::cerr << ini.error() << std::endl;
        exit(-1);
    }
    auto configs = NN::NN_configs(ini);
    std::cout << "----output file is written in: " << configs.result_dir << configs.result_name << "*" << std::endl;

    // set parallel threads in openpm.
    if (configs.use_omp)
    {
        omp_set_num_threads(configs.n_omp_threads);
        std::cout << "----number of threads of openmp: " << configs.n_omp_threads << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    //* main program:
    write_dat(configs);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout.precision(4);
    std::cout << "Duration: " << duration.count() << " seconds\n";
}