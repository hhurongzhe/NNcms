#include "lib_define.hpp"
#include "interaction_all.hpp"

// pack the readable file "txtfname" to binary file "binfname",
// they have the same number of double values.
void pack_kernel_file(const std::string &txtfname, const std::string &binfname, const size_t num)
{
    std::ifstream fp_txt(txtfname);
    if (!fp_txt.is_open())
    {
        std::cerr << "failed to open interaction file: " << txtfname << std::endl;
    }
    std::ofstream fp_bin(binfname, std::ios::binary);
    for (size_t i = 0; i < num; i = i + 1)
    {
        double value;
        fp_txt >> value;
        fp_bin.write(reinterpret_cast<const char *>(&value), sizeof(double));
    }
    fp_txt.close();
    fp_bin.close();
}

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

    // write txt file for this partial-wave channel.
    std::ostringstream oss_txt;
    oss_txt << configs.result_dir << "kernel-" << configs.result_name << "-" << l_final << "-" << l_initial << "-" << s << "-" << j << "-" << tz_name << ".txt";
    auto file_txt_name_this_channel = oss_txt.str();
    std::cout << "writing: " << file_txt_name_this_channel << std::endl;
    std::ofstream fp(file_txt_name_this_channel);
    if (!fp.is_open())
    {
        std::cerr << "failed to open file: " << file_txt_name_this_channel << "!\n";
        exit(-1);
    }
    for (size_t idx_mom_bra = 0; idx_mom_bra < configs.mesh_points_number; idx_mom_bra = idx_mom_bra + 1)
    {
        for (size_t idx_mom_ket = 0; idx_mom_ket < configs.mesh_points_number; idx_mom_ket = idx_mom_ket + 1)
        {
            double p_final = configs.momentum_mesh_points[idx_mom_bra];
            double p_initial = configs.momentum_mesh_points[idx_mom_ket];
            double v_value = interaction_all::potential_chiral(l_final, l_initial, s, j, tz, p_final, p_initial, configs);
            fp << std::fixed << " " << std::scientific << std::setprecision(17) << v_value;
        }
        fp << "\n";
    }
    fp.close();

    // pack binary file from txt file.
    std::ostringstream oss_bin;
    oss_bin << configs.result_dir << "kernel-" << configs.result_name << "-" << l_final << "-" << l_initial << "-" << s << "-" << j << "-" << tz_name << ".bin";
    auto file_bin_name_this_channel = oss_bin.str();
    pack_kernel_file(file_txt_name_this_channel, file_bin_name_this_channel, configs.mesh_points_number * configs.mesh_points_number);
    std::cout << "writing: " << file_bin_name_this_channel << std::endl;
}

// write results in the output file.
void write_dat(const NN::NN_configs &configs)
{
    // write momentum mesh.
    std::ostringstream oss_mom_mesh;
    oss_mom_mesh << configs.result_dir << configs.result_name << "-momentum-mesh.txt";
    auto file_mom_mesh = oss_mom_mesh.str();
    std::ofstream fp_mom_mesh(file_mom_mesh);
    fp_mom_mesh << "# momentum mesh points number;\n";
    fp_mom_mesh << "# momentum mesh points and weights;\n";
    fp_mom_mesh << configs.mesh_points_number << "\n";
    for (size_t i = 0; i < configs.momentum_mesh_points.size(); i = i + 1)
    {
        fp_mom_mesh << std::fixed << std::scientific << std::setprecision(17) << configs.momentum_mesh_points[i] << "\t" << std::setw(17) << configs.momentum_mesh_weights[i] << std::endl;
    }
    fp_mom_mesh.close();

    // write partial-waves.
    std::ostringstream oss_pws;
    oss_pws << configs.result_dir << configs.result_name << "-partial-waves.txt";
    auto file_pws = oss_pws.str();
    std::ofstream fp_pws(file_pws);
    fp_pws << "# partial-waves: l' l s j tz;\n";
    for (size_t i = 0; i < configs.partial_waves.size(); i = i + 1)
    {
        fp_pws << configs.partial_waves[i][0] << " " << configs.partial_waves[i][1] << " " << configs.partial_waves[i][2] << " " << configs.partial_waves[i][3] << " " << configs.partial_waves[i][4] << "\n";
    }
    fp_pws.close();

    // generate channels.
    auto channels = configs.partial_waves;
    for (size_t idx_channel = 0; idx_channel < channels.size(); idx_channel = idx_channel + 1)
    {
        // write matrix elements for each channel.
        auto this_channel = channels[idx_channel];
        write_dat_single_channel(this_channel, configs);
    }
}

int main()
{
    std::cout << "---- running NN-cms...\n\n";

    //---- print current date.
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    char dateStr[100];
    std::strftime(dateStr, sizeof(dateStr), "%Y-%m-%d", std::localtime(&now_c));
    std::cout << "---- current Date: " << dateStr << "\n"
              << std::endl;

    //---- config file initializing.
    auto ini = inifile_system::inifile("inifile-cms.ini");
    if (!ini.good())
    {
        std::cerr << ini.error() << std::endl;
        exit(-1);
    }
    auto configs = NN::NN_configs(ini);
    std::cout << "---- output file is written in: " << configs.result_dir << "kernel-" << configs.result_name << "-ll-l-s-j-tzname.txt & .bin" << std::endl;
    std::cout << "                                " << configs.result_dir << configs.result_name << "-momentum-mesh.txt" << std::endl;
    std::cout << "                                " << configs.result_dir << configs.result_name << "-partial-waves.txt\n"
              << std::endl;

    //---- set parallel threads in openpm, shouldn't too large.
    const size_t thread_number = 16;
    std::cout << "---- number of threads for openmp: " << thread_number << "\n"
              << std::endl;
    omp_set_num_threads(thread_number);

    auto start = std::chrono::high_resolution_clock::now();

    //---- main program:
    write_dat(configs);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout.precision(4);
    std::cout << "\nduration: " << duration.count() << " seconds\n";
}