#ifndef LINEAR_MODE_H
#define LINEAR_MODE_H

#include <vector>


void linearFit(std::vector<double>& x, std::vector<double>& y, double& m, double& q);


void bisezione_lin(std::vector<double>& par, std::vector<double> m, std::vector<double> q);


// Funzione per trovare il parametro "peggiore" -> devo capire qual è il parametro peggiore in funzione del quale migliorare gli altri con 'linear_mode'
int worst_parameter(std::vector<double>& m_lin, std::vector<double>& q_lin);

//funzione effettiva da usare nella main
void linear_mode(std::vector<std::vector<double>>& par_lin, std::vector<double>& m_lin, std::vector<double>& q_lin, bool& errore_lin, bool& ricerca_retta, bool& faster, bool complex);



#endif