#include "automatic_research.h"
#include "covering.h"
#include "interpolating_function.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

extern double chi_quadro_min;       //Chi quadro minimo assoluto


//Ricerca automatica logaritmica
void ricerca_auto(vector<double>& par, vector<double>& par_auto, vector<double>& par_def, int n) {
    //mettere 'par_auto=par' e 'par_auto=par' prima della chiamata della funzione

    static double chi_min_auto = 1e25;
    static vector<double> sec_best_par = par;
    static double min = -1e10;
    static double max = 1e10;
    static double min_ordine = 0.001;
    static bool salita = true;         //Sto scorrendo in "salita" i parametri fino ad arrivare all'ultimo per poi cominciare a tornare indietro per migliorarli uno a uno "all'indietro"
    static int par_size = par.size();

    if (n >= par_size) {        // || (par_auto[n] != 0 && par[n] == 0)
        salita = false;
        return;
    }


    double a = min;
    while (a <= max)
    {
        //cout << "a = " << a << "e il max = " << max << endl;
        //if (a <= max) cout << "ma che cazz..." << endl;
        int k = 0;
        if (n < par_size)
        {
            if (par[n] == 0) {
                par_auto[n] = a;
                //n++;
                //cout << "messo pari ad a e con n = " << n << endl;
            }
            else {          //if (salita)
                while (n + k < par_size && par[n + k] != 0)
                    k++;
                //cout << "Sono a n = " << n << " e k = " << k << endl;
                k--;
            }
        }

        //if (salita)
        ricerca_auto(par, par_auto, par_def, n + 1 + k);

        //Se sono in 'salita', qui sono arrivato all'ultimo parametro e ora comincio a tornare indietro a migliorare quelli precedenti (se sono parametri automatici "a")

        if (par[n] != 0) return;    //Torno indietro se ne trovo uno che non è di quelli automatici "a" e procedo con quello successivo

        double sum_chi = f_chi_quadro(par_auto);

        //Debug
        if (false) {
            cout << "n = " << n << endl;
            for (int k = 0; k < par_size; k++)
            {
                cout << par_auto[k] << "\t";
            }
            cout << endl;
            cout << sum_chi << endl;
        }


        if (sum_chi <= chi_min_auto)
        {
            sec_best_par = par_def;         // Il secondo miglior set di parametri
            par_def = par_auto;             // Il miglior set di parametri
            chi_min_auto = sum_chi;
            /*
            if (min_ordine > 1e-5 && (a >= -min_ordine - 1e-9 && a <= -min_ordine + 1e-9))
                min_ordine /= 10;
            */
            //cout << " --> yep" << endl;
        }

        // Devo andare da -1e10 a 1e10 togliendo [-1e-2;+1e-2]
        (a < 0) ? a /= 10 : a *= 10;
        if (a >= -min_ordine / 10 - 1e-9 && a <= -min_ordine / 10 + 1e-9)
            a = min_ordine;
        //cout << "a = " << a << endl;
    }
    //cout << "n = " << n << "e sono uscito dal while" << endl;

    //Sono uscito dal while e verifico se sono proprio alla fine
    auto it = find(par.begin(), par.end(), 0);      // restituisce un puntatore al primo elemento uguale a 0 (parametro automatico)
    int index = distance(par.begin(), it);          // restituisce l'indice di tale elemento

    if (n == index)     //Sono alla fine di tutta l'esecuzione della funzione
    {
        par = par_def;              // Alla fine 'par'      sono i migliori parametri
        par_auto = sec_best_par;    // Alla fine 'par_auto' sono i secondi migliori parametri
    }
}


void parametri_auto(vector<double>& par, bool output) {

    //Ricerca automatica logaritmica
    vector<double> par_auto = par;
    vector<double> par_def = par;
    ricerca_auto(par, par_auto, par_def, 0);

    chi_quadro_min = f_chi_quadro(par);
    if (output)
    {
        cout << "Parametri da ricerca automatica:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        cout << endl << "Chi quadro = " << chi_quadro_min << endl;
    }


    //Cerco di capire se è meglio 'par' o 'par_auto' con algoritmo di bisezione ('par_auto' è il secondo migliore)
    vector<double> passo_i1;
    vector<double> passo_i2;

    vector<double> par_i1(par.size(), 0);
    vector<double> par_i2(par.size(), 0);

    algoritmo_bisezione(par, par_i1, passo_i1, 0);
    algoritmo_bisezione(par_auto, par_i2, passo_i1, 0);

    par = (f_chi_quadro(par_i1) < f_chi_quadro(par_i2)) ? par_i1 : par_i2;

    chi_quadro_min = f_chi_quadro(par);
    if (output)
    {
        cout << "Parametri automatici iniziali migliorati:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        cout << endl << "Chi quadro = " << chi_quadro_min << endl;
    }
}