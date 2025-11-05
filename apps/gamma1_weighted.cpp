// gamma1 weighted (equation 29 implementation)
#include "orbit.h"
#include "gamma.h"
#include "io.h"
#include "precision.h"

#include <iostream>
#include <fstream>
#include <iomanip>

int main() {
    const std::string input_file = "input/gamma1_weighted.input";
    const std::string output_file = "output/gamma1_weighted.output";
    
    auto params = parse_input(input_file);

    OrbitParams p = {
        params["Mp"], params["Mc"], params["Porb"], params["i"],
        params["e"], 0.0L, params["omega"]
    };

    real T_obs = params["T"];
    int m = static_cast<int>(params["m"]);
    real Pp_s = params["Pspin"];

    std::cout << "Starting orbit-averaged gamma1 computation..." << std::endl;

    real sum_num = 0.0L;
    real sum_den = 0.0L;
	const int f_steps = 36; // 10 degree step-size
    // const int f_steps = 360; // 1 degree step-size (use this for more accurate values)

    // #pragma omp parallel for reduction(+:sum_num, sum_den)
    for (int i = 0; i < f_steps; ++i) {
        real f0_deg = 360.0L * ((real)i) / ((real)f_steps); 

        
        real gamma_f = gamma1(p, f0_deg, T_obs, m, Pp_s);
        real w_f = compute_w(f0_deg, p.e);

        sum_num += gamma_f * w_f;
        sum_den += w_f;
    }

    real gamma_avg = sum_num / sum_den;
     
    std::cout << "Computation finished." << std::endl;
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "  Weighted gamma1 = " << gamma_avg << std::endl;

    // write to output file
    std::ofstream outfile(output_file);
    outfile.setf(std::ios::fixed);
    outfile << std::setprecision(8);
    outfile << "# Weighted (orbit-averaged) gamma1 calculation result" << std::endl;
    outfile << "gamma1_weighted " << gamma_avg << std::endl;
    outfile.close();

    std::cout << "Result written to " << output_file << std::endl;

    return 0;
}

