#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>

bool isDataLine(const std::string& line);

double findMinIgnoringZero(const std::vector<double>& vec);


int readFile(const std::string filePath, std::vector<double>& x, std::vector<double>& sigma_x, std::vector<double>& y, std::vector<double>& sigma_y);

void writeFile(const std::string filePath, std::vector<double>& x, std::vector<double>& sigma_x, std::vector<double>& y, std::vector<double>& sigma_y, std::vector<double> par_best, bool approx, std::string function);


#endif