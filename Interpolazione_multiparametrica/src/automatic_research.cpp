#include "automatic_research.h"

#include "interpolator.h"
#include "gradient_descent_algorithm.h"
#include "BranchingCover.h"
//#include "ThreadPool.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

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

    //Cerco di capire qual è il miglior 'par' con 'BranchingCover::doBetter'
    int n = (par_improved.size() < 5) ? par_improved.size() : 5;
    for (int i = 0; i < n; i++)
    {
        std::vector<double> par_prov = par_improved[par_improved.size() - 1 - i];
        double chi_min = (i == 0) ? chi_quadro_min : i_generator->fChiQuadro(par_prov);
        if (chi_min > chi_quadro_min * 100)     // se il rispettivo chi_quadro è più di 100 volte più grande di quello minimo delle altre scelte non li considero più da lì in poi
            return;

        // Metodo di ramificazione per migliorare i parametri
        BranchingCover b_generator(par_prov);
        b_generator.compute(par_prov);

        chi_min = i_generator->fChiQuadro(par_prov);
        if (chi_min < chi_quadro_min) {
            par = par_prov;
            chi_quadro_min = chi_min;
            if (output) print("Auto-search parameters iteration " + std::to_string(i), par);
        }
        else
            break;

    }

    std::vector<double> par_prov = par;
    BranchingCover b_generator2(par_prov);
    b_generator2.compute(par_prov);
    double chi_min = i_generator->fChiQuadro(par_prov);
    if (chi_min < chi_quadro_min) {
        par = par_prov;
        chi_quadro_min = chi_min;
        if (output) print("Auto-search parameters", par);
    }
    

}


void AutomaticResearch::endJob(std::vector<double>& par_best) {
    par_best = par;
    if (output) print("Improved auto-search parameters:", par_best);
}


void AutomaticResearch::print(const std::string name, std::vector<double> par) {
    std::cout << name << std::endl;
    for (int k = 0; k < par.size(); k++)
        std::cout << par[k] << "\t";

    std::cout << std::endl << "Chi-squared = " << i_generator->fChiQuadro(par) << std::endl;
}



AutomaticResearch::~AutomaticResearch() {
    delete i_generator;
}