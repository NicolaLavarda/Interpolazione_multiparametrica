#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <string>

using namespace std;

//double funzione_interpolante(vector<double> x, vector<double> par, int i);


void ricerca_auto(vector<double>& par, vector<double>& par_auto, int n);

void ricerca_auto(vector<double>& par, vector<double>& par_auto, vector<double>& par_def, int n);

void algoritmo_bisezione(vector<double> par, vector<double>& par_def, const vector<double> passo, int n);

void ricoprimento(vector<double>& par, vector<double>& par_def, const vector<double> passo, int livello, int dimensione, bool is_on_surface);

double chi_quadro_piu_uno(int num_parametri);

double p_value(double chi_observato, int GDL);

double chi_piu_uno_par_n(vector<double> parametri, int numero_parametro, int chi_piu_uno_val);

void risultato(double valore, double errore, string nome, bool arrotondamento);

void risultato1(vector<double> par_best, bool approx_bool);

void risultato2(vector<double> par_best, bool approx_bool);

#endif