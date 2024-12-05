#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>

int readFile(const std::string filePath, std::vector<double>& x, std::vector<double>& sigma_x, std::vector<double>& y, std::vector<double>& sigma_y);

//void Popolamento_vettori(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sigma_y, std::string file_dati);


#endif