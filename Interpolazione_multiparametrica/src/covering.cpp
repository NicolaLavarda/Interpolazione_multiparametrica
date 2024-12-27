#include "covering.h"
#include "interpolating_function.h"
#include "bisection_algorithm.h"
#include "matrix.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>  // per std::max_element

using namespace std;


extern double chi_quadro_min;                                    //Chi quadro minimo assoluto
extern vector<double> chi_quadro;                                //Vettore chi quadri minimi
extern vector<double> par_best;
extern vector<int> livelli;
extern vector<vector<double>> par_matrix;
extern int cicle_programms;                                                //numero di cicli di riduzione di 'spostamento'


// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void ricoprimento(vector<double>& par, vector<double>& par_def, const vector<double> passo, int livello, int dimensione, bool is_on_surface, bool output) {
    vector<double> par_prov = par_matrix[cicle_programms - 1]; //non va bene 'static vector<double> par_prov = par;' perché altrimenti mi rimane anche per le "riesecuzioni del programma successive". Importante che sia esattamente così: se passassi in qualsiasi modo il vettore alla funzione fidati che non funziona perché poi cambia, fidati.

    if (dimensione == par.size()) {
        if (is_on_surface) {

            //Algoritmo di bisezione --------------------------------------------------------------------
            
            double sum_chi_p = 0;
            vector<double> par_bis = par;                       // 'par' ha le coordinate giuste del centro del cubo n-dimensionale in analisi, quindi dev'essere l'argomento della funzione, ma non può essere modificato (cosa che invece farebbe la funzione), quindi gli passo una copia
            bisection_algorithm(par_bis, passo, sum_chi_p);
            par_def = par_bis;      // 'par_def' sono ora quelli migliorati con l'algoritmo di bisezione

            //cout << par_def[0] << "\t" << par_def[1] << "\t" << sum_chi_p << endl;

            if (sum_chi_p < chi_quadro_min) {
                /*
                    string output_testo = "Optimized parameters: \n";
                    for (int p = 0; p < par.size(); p++) {
                        output_testo += "par" + to_string(p) + ":\t" + to_string(par_def[p]) + "\n";
                    }
                    output_testo += "chi_quadro = " + to_string(sum_chi_p) + "\n";
                    cout << output_testo << endl;           //output valori
                */
                if (output)
                    cout << "\t" << sum_chi_p << " ";

                chi_quadro_min = sum_chi_p;
                par_best = par_def;
                chi_quadro.push_back(sum_chi_p);          //Lo metto nella main in modo da aggiungerlo solo alla fine di un livello di ricoprimento (in uno stesso livello può essere che rimanga miniore del "per mille")
                livelli[livello]++;                       //Metto un valore diverso da 0 se in un livello riesco a migliorare il chi quadro
            }
            //-------------------------------------------------------------------------------------------
        }
        return;
    }

    // Generiamo tutti gli spostamenti per questa dimensione usando il passo specifico (double)
    for (double delta = -livello * passo[dimensione]; delta <= (livello * passo[dimensione] + 1e-9); delta += passo[dimensione]) {
        int m = 0;
        for (int t = 0; t < par.size(); t++)
        {
            if ((fabs(par[t] - (par_prov[t] - livello * passo[t])) < 1e-9) || (fabs(par[t] - (par_prov[t] + livello * passo[t])) < 1e-9)) m++;
        }
        bool surface_condition = m;
        //bool surface_condition = (fabs(delta - livello * passo[dimensione]) < 1e-9);  // Precisione double circa pari a 1e-9
        par[dimensione] = par_prov[dimensione] + delta;  // Aggiorna la coordinata corrente
        ricoprimento(par, par_def, passo, livello, dimensione + 1, is_on_surface || surface_condition, output);
    }
}

