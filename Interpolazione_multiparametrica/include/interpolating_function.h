#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>
#include <string>

using namespace std;


void setup_expression(string expression_str);

double funzione_interpolante(vector<double>& x, vector<double>& par, int i);

double funzione_interpolante1(vector<double> x, vector<double> par, int i);		//non serve a nulla, è per sicurezza se non funzionasse 'funzione_interpolante'

double f_chi_quadro(vector<double> par);


#endif
