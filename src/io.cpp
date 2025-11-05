#include "precision.h"
#include "io.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>

std::map<std::string, real> parse_input(const std::string& filename) {
    std::map<std::string, real> params;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        std::exit(1);
    }

    std::string line, key;
    real value;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // skip comments
        std::stringstream ss(line);
        if (ss >> key >> value) {
            params[key] = value;
        }
    }
    return params;
}
