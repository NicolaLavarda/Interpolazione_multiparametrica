#ifndef FILE_H
#define FILE_H

#include "input.h"

#include <vector>
#include <string>


int readFile(const std::string filePath, std::vector<double>& x, std::vector<double>& sigma_x, std::vector<double>& y, std::vector<double>& sigma_y);

void writeFile(input::Data data, input::Interpolation interpolation, bool approx);


#endif