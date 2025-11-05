#ifndef IO_H
#define IO_H

#include "precision.h"
#include <string>
#include <map>

// Function to parse a simple "key value" input file
std::map<std::string, real> parse_input(const std::string& filename);

#endif // IO_H