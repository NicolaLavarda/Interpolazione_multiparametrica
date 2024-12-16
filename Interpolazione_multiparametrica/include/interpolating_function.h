#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>
#include <string>


void setup_expression(std::string expression_str);

double dfdx(int i);

double f_chi_quadro(std::vector<double> par);

double funzione_interpolante(std::vector<double> par, double x);

double x_function(std::vector<double> par, int i);

#endif
