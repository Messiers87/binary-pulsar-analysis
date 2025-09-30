// gamma_weightedLooped 
#include "../src/orbit.h"
#include "../src/gamma.h"
#include "../src/io.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

int main() {
    const std::string input_file = "input/gamma1_weighted_loopPorbPspin.input";
    const std::string output_file = "output/gamma1_weighted_loopPorbPspin.output";
    
    auto params = parse_input(input_file);
    
    OrbitParams base_p = {
        params["Mp"], params["Mc"], 0.0L,
        params["i"], params["e"], 0.0L, params["omega"]
    };

    long double T_obs = params["T"];
    int m = static_cast<int>(params["m"]);

    std::vector<long double> Po_days_vec = { 0.1L, 0.2L, 0.5L, 1.0L, 10.0L, 100.0L};
    std::vector<long double> Pp_s_vec = { 0.001L, 0.005L, 0.01L, 0.05L, 0.1L, 1.0L };

    std::ofstream outfile(output_file);
    outfile.setf(std::ios::fixed);
    outfile << std::setprecision(8);
    outfile << "# Porb_day      Pspin_s       gamma1_weighted" << std::endl;

    for (auto Po_d : Po_days_vec) {
        for (auto Pp_s : Pp_s_vec) {
            OrbitParams current_p = base_p;
            current_p.Po_day = Po_d;

            std::cout << "Processing: Porb = " << Po_d << " days, Pspin = " << Pp_s << " s" << std::endl;
            
            long double sum_num = 0.0L;
            long double sum_den = 0.0L;
            const int f_steps = 180;

            const long double w_max_denom = powl(1.0L + current_p.e * cosl(deg2rad(180.0L)), -2.0L);

            #pragma omp parallel for reduction(+:sum_num, sum_den)
            for (int i = 0; i < f_steps; ++i) {
                long double f0_deg = 360.0L * ((long double)i) / ((long double)f_steps);
                
                long double gamma_f = gamma1(current_p, f0_deg, T_obs, m, Pp_s);
                long double w_num = powl(1.0L + current_p.e * cosl(deg2rad(f0_deg)), -2.0L);
                long double w_f = w_num / w_max_denom;

                sum_num += gamma_f * w_f;
                sum_den += w_f;
            }
            long double gamma_avg = sum_num / sum_den;

            outfile << std::setw(14) << Po_d << " " << std::setw(14) << Pp_s << " " << std::setw(14) << gamma_avg << std::endl;
            std::cout << "  -> Weighted gamma1 = " << gamma_avg << std::endl;
        }
    }

    outfile.close();
    std::cout << "\nAll computations finished. Results written to " << output_file << std::endl;

    return 0;
}