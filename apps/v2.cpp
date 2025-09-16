// gamma1 weighted (equation 29 implementation)
// produce csv data with columns gamma1avg and w
#include "../src/orbit.h"
#include "../src/gamma.h"
#include "../src/io.h"

#include <iostream>
#include <fstream>
#include <iomanip>

// int main() {
//     auto [params, run] = read_input("../input/gamma1_weighted.input");

//     long double T_obs = 500.0L;
//     int m = 4;
//     long double Pp_s = 0.01L;

//     long double sum_num = 0.0L;
//     long double sum_den = 0.0L;

//     for (int fdeg = 0; fdeg < 360; fdeg += 2) {
//         long double gamma_f = gamma1(params, fdeg, T_obs, m, Pp_s);
//         long double w = compute_w(fdeg, params.e);
//         sum_num += gamma_f * w;
//         sum_den += w;
//     }

//     long double gamma_avg = sum_num / sum_den;

//     write_output_single("../output/gamma1_weighted.output", gamma_avg );
//     std::cout << "Weighted gamma1=" << gamma_avg << " written to gamma1_weighted.output\n";

//     return 0;
// }

int main() {
    const std::string input_file = "input/gamma1_weighted.input";
    const std::string output_file = "output/gamma1_weighted.output";
    
    auto params = parse_input(input_file);

    OrbitParams p = {
        params["Mp"], params["Mc"], params["Porb"], params["i"],
        params["e"], 0.0L, params["omega"]
    };

    long double T_obs = params["T"];
    int m = static_cast<int>(params["m"]);
    long double Pp_s = params["Pspin"];

    std::cout << "Starting orbit-averaged gamma1 computation..." << std::endl;

    long double sum_num = 0.0L;
    long double sum_den = 0.0L;
    const int f_steps = 180;

    #pragma omp parallel for reduction(+:sum_num, sum_den)
    for (int i = 0; i < f_steps; ++i) {
        long double f0_deg = 360.0L * ((long double)i) / ((long double)f_steps);
        
        long double gamma_f = gamma1(p, f0_deg, T_obs, m, Pp_s);
        long double w_f = powl(1.0L + p.e * cosl(deg2rad(f0_deg)), -2.0L);

        sum_num += gamma_f * w_f;
        sum_den += w_f;
    }

    long double gamma_avg = sum_num / sum_den;
    
    std::cout << "Computation finished." << std::endl;
    std::cout << "  Weighted gamma1 = " << gamma_avg << std::endl;
    
    std::ofstream outfile(output_file);
    outfile.setf(std::ios::fixed);
    outfile << std::setprecision(8);
    outfile << "# Weighted (orbit-averaged) gamma1 calculation result" << std::endl;
    outfile << "gamma1_weighted " << gamma_avg << std::endl;
    outfile.close();

    std::cout << "Result written to " << output_file << std::endl;

    return 0;
}
