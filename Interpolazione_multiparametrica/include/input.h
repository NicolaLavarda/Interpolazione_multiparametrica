#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <string>
#include <map>

void input(int argc, char* argv[], std::vector<double>& par,
		   int& num_a, std::map<std::string, bool>& options,
		   std::vector<double>& x, std::vector<double>& sigma_x,
		   std::vector<double>& y, std::vector<double>& sigma_y);

bool improved(const std::string filePath, std::vector<double> par_best);

double GetChiSquared(const std::string filePath);

#endif