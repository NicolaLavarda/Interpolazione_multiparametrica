#ifndef CHI_SQUARE_H
#define CHI_SQUARE_H

#include <vector>


double chi_quadro_piu_uno(int num_parametri);

double p_value(double chi_observato, int GDL);

double chi_piu_uno_par_n(std::vector<double> parametri, int numero_parametro, int chi_piu_uno_val);


#endif