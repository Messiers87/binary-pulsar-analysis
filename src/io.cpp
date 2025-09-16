#include "io.h"
#include <fstream>
#include <iostream>
#include <sstream>

std::map<std::string, long double> parse_input(const std::string& filename) {
    std::map<std::string, long double> params;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        exit(1);
    }
    std::string line, key;
    long double value;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        std::stringstream ss(line);
        if (ss >> key >> value) {
            params[key] = value;
        }
    }
    return params;
}