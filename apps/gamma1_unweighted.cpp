// gamma1_unweighted
#include "../src/orbit.h"
#include "../src/gamma.h"
#include "../src/io.h"

#include <iostream>
#include <fstream>
#include <iomanip>

constexpr long double PI = 3.141592653589793238462643383279502884L;

int main() {
    const std::string input_file = "input/gamma1_unweighted.input";
    const std::string output_file = "output/gamma1_unweighted.output";

    auto params = parse_input(input_file);

    OrbitParams p = {
        params["Mp"], params["Mc"], params["Porb"], params["i"],
        params["e"], 0.0L, params["omega"]
    };

    long double f0_deg = params["f0"];
    long double T_obs = params["T"];
    int m = static_cast<int>(params["m"]);
    long double Pp_s = params["Pspin"];

    std::cout << "Starting unweighted gamma1 computation for f0 = " << f0_deg << " deg..." << std::endl;

    long double gamma = gamma1(p, f0_deg, T_obs, m, Pp_s);
    long double num = powl(1.0L + p.e * cosl(deg2rad(f0_deg)), -2.0L);
    long double denom = powl(1.0L + p.e * cosl(PI), -2.0L);
	long double w = num / denom;
    return w;

    std::cout << "Computation finished." << std::endl;
    std::cout << "  gamma1 = " << gamma << std::endl;
    std::cout << "  w      = " << w << std::endl;

    std::ofstream outfile(output_file);
    outfile.setf(std::ios::fixed);
    outfile << std::setprecision(8);
    outfile << "# Unweighted gamma1 calculation results" << std::endl;
    outfile << "gamma1 " << gamma << std::endl;
    outfile << "w      " << w << std::endl;
    outfile.close();

    std::cout << "Results written to " << output_file << std::endl;

    return 0;
}