#ifndef COVERING_H
#define COVERING_H

#include <vector>


void algoritmo_bisezione(std::vector<double> par, std::vector<double>& par_def, const std::vector<double> passo, int n);

void ricoprimento(std::vector<double>& par, std::vector<double>& par_def, const std::vector<double> passo, int livello, int dimensione, bool is_on_surface);


std::vector<double> grad_f_chi_quadro(const std::vector<double> par);

void discesa_gradiente(std::vector<double>& par, double sensibility);


#endif