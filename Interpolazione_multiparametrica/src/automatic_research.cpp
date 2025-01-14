#include "automatic_research.h"

#include "interpolator.h"
#include "gradient_descent_algorithm.h"
//#include "ThreadPool.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

AutomaticResearch::AutomaticResearch(std::vector<double> par, bool output):
    par(par), output(output)
{
    auto it = find(par.begin(), par.end(), 0);      // restituisce un puntatore al primo elemento uguale a 0 (parametro automatico)
    first = distance(par.begin(), it);          // restituisce l'indice di tale elemento
}


void AutomaticResearch::compute(std::vector<double>& par_auto, int n) {
    static double chi_min_auto = 1e25;
    static double min = -1e10;
    static double max = 1e10;
    static double min_ordine = 0.001;       //dev'essere positivo (>0)
    static bool salita = true;         //Sto scorrendo in "salita" i parametri fino ad arrivare all'ultimo per poi cominciare a tornare indietro per migliorarli uno a uno "all'indietro"
    static int par_size = par_auto.size();

    if (n >= par_size) {
        salita = false;
        return;
    }

    double a = min;
    while (a <= max)
    {
        int k = 0;
        if (n < par_size)
        {
            if (par[n] == 0) {
                par_auto[n] = a;
            }
            else {          //if (salita)
                while (n + k < par_size && par[n + k] != 0)
                    k++;
                k--;
            }
        }

        //if (salita)
        compute(par_auto, n + 1 + k);

        //Se sono in 'salita', qui sono arrivato all'ultimo parametro e ora comincio a tornare indietro a migliorare quelli precedenti (se sono parametri automatici "a")

        if (par[n] != 0) return;    //Torno indietro se ne trovo uno che non è di quelli automatici "a" e procedo con quello successivo

        double sum_chi = i_generator.fChiQuadro(par_auto);

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
    
    double chi_quadro_min = i_generator.fChiQuadro(par);
    if (output)
    {
        cout << "Parametri da ricerca automatica:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        cout << endl << "Chi quadro = " << chi_quadro_min << endl;
    }

    // Creo un pool di multi-thread
    //ThreadPool pool;

    // Aggiungo task al pool
    int num_it = (par_improved.size() < 10) ? par_improved.size() : 10;
    for (int i = 0; i < num_it; i++)
    {
        //analysis.emplace_back(pool.enqueue(discesa_gradiente, i, 10));
    }

    // Recupero i risultati
    for (int i = 0; i < num_it; i++) {
        //std::cout << "Risultato: " << analysis[i].get() << std::endl;
    }

    //Cerco di capire se è meglio 'par' o 'par_i2' con "discesa_gradiente" ('par_i2' è il secondo migliore)
    std::vector<double> par_i1 = par;                                   // le due alternative
    std::vector<double> par_i2 = par_improved[par_improved.size()-2];   //

    double sensibility0 = 1;    // variazione iniziale del 100%
    gradient_descent_algorithm(par_i1, chi_quadro_min, sensibility0);    //Miglioro 'par_i1' con metodo "discesa_gradiente"

    double sensibility1 = 1;    // variazione iniziale del 100%
    gradient_descent_algorithm(par_i2, chi_quadro_min, sensibility1);    //Miglioro 'par_i2' con metodo "discesa_gradiente"

    par = (i_generator.fChiQuadro(par_i1) < i_generator.fChiQuadro(par_i2)) ? par_i1 : par_i2;      //Controllo se è meglio 'par_i1' o 'par_i2'
    chi_quadro_min = i_generator.fChiQuadro(par);

    if (output)
    {
        cout << "Parametri automatici iniziali migliorati:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        cout << endl << "Chi quadro = " << chi_quadro_min << endl;
    }
}


void AutomaticResearch::endJob(std::vector<double>& par_best) {
    par_best = par;
}