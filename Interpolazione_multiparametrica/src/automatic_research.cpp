#include "automatic_research.h"

#include "interpolator.h"
#include "gradient_descent_algorithm.h"
//#include "ThreadPool.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

AutomaticResearch::AutomaticResearch(std::vector<double> par, bool output):
    par(par), output(output), par_size(par.size())
{
    auto it = find(par.begin(), par.end(), 0);      // restituisce un puntatore al primo elemento uguale a 0 (parametro automatico)
    first = distance(par.begin(), it);          // restituisce l'indice di tale elemento
}


void AutomaticResearch::compute(std::vector<double>& par_auto, int n) {    

    if (n >= par_size) return;

    double a = min;
    while (a <= max)
    {
        int k = 0;
        if (n < par_size)
        {
            if (par[n] == 0) {
                par_auto[n] = a;
            }
            else {
                while (n + k < par_size && par[n + k] != 0)
                    k++;
                k--;
            }
        }

        compute(par_auto, n + 1 + k);

        //Se sono in "salita", qui sono arrivato all'ultimo parametro e ora comincio a tornare indietro a migliorare quelli precedenti (se sono parametri automatici "a")

        if (par[n] != 0) return;    //Torno indietro se ne trovo uno che non è di quelli automatici "a" e procedo con quello successivo

        double sum_chi = i_generator->fChiQuadro(par_auto);

        if (sum_chi <= chi_min_auto)
        {
            par_improved.push_back(par_auto);         // parametri di volta in volta migliorati
            chi_min_auto = sum_chi;
        }

        // Devo andare da -1e10 a 1e10 togliendo [-1e-2;+1e-2]
        (a < 0) ? a /= 10 : a *= 10;
        if (a >= -min_ordine / 10 - 1e-9 && a <= -min_ordine / 10 + 1e-9)
            a = min_ordine;
    }

    //Sono uscito dal while e verifico se sono proprio alla fine
    if (n == first)
    {
        //Sono alla fine di tutta l'esecuzione della funzione
    }
}


void AutomaticResearch::beginJob() {
    //Ricerca automatica logaritmica
    std::vector<double> par_auto = par;
    compute(par_auto, 0);
    par = par_improved.back();
    
    double chi_quadro_min = i_generator->fChiQuadro(par);

    if (output) print("Auto-search parameters:", par);

    //Cerco di capire qual è il miglior 'par' con "discesa_gradiente"
    int n = (par_improved.size() < 5) ? par_improved.size() : 5;
    double sensibility = 0.1;
    for (size_t i = 0; i < n; i++)
    {
        par = par_improved[par_improved.size() - 1 - i];
        gradient_descent_algorithm(par, chi_quadro_min, sensibility);
    }

    if (output) print("Improved auto-search parameters:", par);

}


void AutomaticResearch::endJob(std::vector<double>& par_best) {
    par_best = par;
}


void AutomaticResearch::print(const std::string name, std::vector<double> par) {
    cout << name << endl;
    for (int k = 0; k < par.size(); k++)
        cout << par[k] << "\t";

    cout << endl << "Chi-squared = " << i_generator->fChiQuadro(par) << endl;
}




void AutomaticResearch::doBetter() {

    // Mi assicuro di partire da un minimo locale
    double sensibility = 0.1;
    double chi_quadro_min = 0.0;
    gradient_descent_algorithm(par, chi_quadro_min, sensibility);
}


AutomaticResearch::~AutomaticResearch() {
    delete i_generator;
}