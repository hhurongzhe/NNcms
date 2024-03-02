#include "lib_define.hpp"
#include "interaction_all.hpp"

// build partial-wave channels before writing matrix elements,
// each element is {l',l,s,j,t,tz}.
std::vector<std::vector<int>> build_partial_wave_channels(const NN::NN_configs &configs)
{
    std::vector<std::vector<int>> temp;
    for (int tz = -1; tz <= 1; tz = tz + 1)
    {
        // j=0 channels: 1S0 and 3P0.
        std::vector<int> channel_1s0 = {0, 0, 0, 0, 1, tz};
        temp.push_back(channel_1s0);
        std::vector<int> channel_3p0 = {1, 1, 1, 0, 1, tz};
        temp.push_back(channel_3p0);
        // j>=1 channels.
        for (int j_temp = 1; j_temp <= configs.J_max; j_temp = j_temp + 1)
        {
            for (int s = 0; s <= 1; s = s + 1)
            {
                for (int l_bra = j_temp - s; l_bra <= j_temp + s; l_bra = l_bra + 1)
                {
                    for (int l_ket = j_temp - s; l_ket <= j_temp + s; l_ket = l_ket + 1)
                    {
                        if (abs(l_bra - l_ket) % 2 != 0) // check parity
                        {
                            continue;
                        }
                        int t = 1 - ((l_bra + s) % 2); // determine t by (-1)^(l+s+t)=odd.
                        if (t == 0 && tz != 0)         // check isospin
                        {
                            continue;
                        }
                        std::vector<int> this_channel = {l_bra, l_ket, s, j_temp, t, tz};
                        temp.push_back(this_channel);
                    }
                }
            }
        }
    }

    return temp;
}

void write_dat_each(std::vector<int> this_channel, const NN::NN_configs &configs)
{
    int l_final, l_initial, s, j, t, tz;
    l_final = this_channel[0];
    l_initial = this_channel[1];
    s = this_channel[2];
    j = this_channel[3];
    t = this_channel[4];
    tz = this_channel[5];
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
    std::ostringstream oss;
    oss << configs.result_dir << configs.result_name << "-" << l_final << "-" << l_initial << "-" << s << "-" << j << "-" << t << "-" << tz_name << ".dat";
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

// write results in the output file.
void write_dat(const NN::NN_configs &configs)
{
    std::ofstream fp(configs.result_file());

    fp << "# momentum mesh points number: " << configs.mesh_points_number << "\n";
    fp << "! momentum mesh points and weights:\n";
    for (int i = 0; i < configs.momentum_mesh_points.size(); i = i + 1)
    {
        fp << std::fixed << "\t" << std::scientific << std::setprecision(17) << configs.momentum_mesh_points[i] << "\t" << std::setw(17) << configs.momentum_mesh_weights[i] << std::endl;
    }

    // generate channels.
    auto channels = build_partial_wave_channels(configs);
    for (int idx_channel = 0; idx_channel < channels.size(); idx_channel = idx_channel + 1)
    {
        // write matrix elements for each channel.
        auto this_channel = channels[idx_channel];
        write_dat_each(this_channel, configs);
    }
}

int main()
{
    // config file initializing.
    auto ini = inifile_system::inifile("inifile.ini");
    if (!ini.good())
    {
        std::cerr << ini.error() << std::endl;
        exit(-1);
    }
    auto configs = NN::NN_configs(ini);
    std::cout << "----output file is written in: " << configs.result_file() << std::endl;

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