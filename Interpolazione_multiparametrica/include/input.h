#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <string>
#include <map>

bool isNumber(const std::string& str);

void input(int argc, char* argv[], std::vector<double>& par, int& num_a, std::map<std::string, bool>& options);

#endif