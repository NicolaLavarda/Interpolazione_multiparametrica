#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>
#include <string>
#include "exprtk.hpp"

using namespace std;

// Dichiara le variabili globali con extern per l'uso in altri file sorgente
extern exprtk::symbol_table<double> symbol_table;
extern exprtk::expression<double> expression;
extern exprtk::parser<double> parser;

extern double* px;
extern double* pa;
extern double* pb;
extern double* pc;

void setup_expression(const std::string& expression_str);
//void setup_expression(string expression_str, vector<double> x, vector<double> par, int i);

//defined in main.cpp
double funzione_interpolante(std::vector<double>& x, std::vector<double>& par, int i);
//double funzione_interpolante(vector<double> x, vector<double> par, int i);







#endif