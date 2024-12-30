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


// costruttore di base
covering::covering(vector<double>& par_best, double chi_quadro_min, double cicle_programms, bool output, bool faster) :
    par_best(par_best), chi_quadro_min(chi_quadro_min), cicle_programms(cicle_programms), output(output), faster(faster)
{

    par_prov = par_best;
    par = par_best;
    chi_quadro_ciclo_prec = chi_quadro_min;

    livello = 0;                         // Livello iniziale

    spostamento = 0.1;       // es. spostamento=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
    if (cicle_programms > 2) spostamento /= cicle_programms * 2;        //aumentando fino ad un '*2' ('spostamento /= cicle_programms * 2') il tempo d'esecuzione è sempre migliorato, ma da valutare bene come ottimizzare il parametro '*2'

    for (int i = 0; i < par.size(); i++)
        passo.push_back(spostamento * ((fabs(par[i]) > 0.5) ? fabs(par[i]) : 2));        //spostamento * 2 se miniore di 0.5 in modo che il rispettivo 'passo' sia pari a 1 (che nella funzione diventerà un range semilargo 0.5)

    //par_matrix.push_back(par_best);     //Faccio arrivare i par_best di ogni intero ciclo di programma alla funzione 'ricoprimento' come utimo vettore del vettore di vettori par_matrix

    //par_matrix_ext = par_matrix;

}

void covering::status(int cicle_val, int k_val) {

    cicle = cicle_val;
    k = k_val;

    if (!output) { cout << "\rNumber of program iterations: " << cicle_programms << flush; }
    else          { cout << endl << "Livello " << cicle << "." << livello << ":";           }
}

// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void covering::ricoprimento(vector<double>& par, vector<double>& par_best, int dimensione, bool is_on_surface) {
    //vector<double> par_prov = par_matrix[cicle_programms - 1]; //non va bene 'static vector<double> par_prov = par;' perché altrimenti mi rimane anche per le "riesecuzioni del programma successive". Importante che sia esattamente così: se passassi in qualsiasi modo il vettore alla funzione fidati che non funziona perché poi cambia, fidati.

    if (dimensione == par.size()) {
        if (is_on_surface) {

            //Algoritmo di bisezione --------------------------------------------------------------------
            
            double sum_chi_p = 0;
            vector<double> par_bis = par;                       // 'par' ha le coordinate giuste dell'attuale centro del cubo n-dimensionale in analisi, quindi dev'essere l'argomento della funzione, ma non può essere modificato (cosa che invece farebbe la funzione), quindi gli passo una copia
            bisection_algorithm(par_bis, passo, sum_chi_p);
            par_best = par_bis;      // 'par_best' sono ora quelli migliorati con l'algoritmo di bisezione

            //cout << par_best[0] << "\t" << par_best[1] << "\t" << sum_chi_p << endl;

            if (sum_chi_p < chi_quadro_min) {

                if (output)
                    cout << "\t" << sum_chi_p << " ";

                chi_quadro_min = sum_chi_p;
                par_best = par_best;
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
        ricoprimento(par, par_best, dimensione + 1, is_on_surface || surface_condition);
    }
}


void covering::next() {

    //Passaggio al ciclo successivo, altrimenti semplicemente passaggio al livello successivo
    if ((livello > 1) && (livelli[livello] == 0 && (livelli[livello - 1] == 0 || (faster || livello > 4))))      //Per livelli 2, 3 e 4 è necessario che almeno gli ultimi due livelli non migliorino il chi quadro per passare al ciclo successivo [...] if (faster || livello > 4) //dal livello 5 in poi basta che un solo livello non migliori il chi quadro per passare al ciclo successivo (oppure se è attivo 'faster')
    {
        livello = 0;
        if (cicle > 2)               //Per il secondo ciclo rimane 'spostamento' uguale al primo ciclo
            spostamento /= 4;
        par = par_best;

        livelli.clear();
        livelli.assign(100, 0);     //livelli viene ripristinato

        for (int i = 0; i < par.size(); i++)
            passo[i] = (fabs(par[i]) > 0.5) ? fabs(par[i] * spostamento) : fabs(2 * spostamento);   //Riduco 'spostamento' per i livelli successivi (all'interno dello stesso 'cicle_programme')
        cicle++;
    }
    else
        ++livello;      // Espansione del ricoprimento ('spostamento' ovviamente rimane uguale, ricopro con cubi n-dimensionali di uguali dimensioni in uno stesso ciclo di ricoprimento)

}

bool covering::exit(bool ricerca_retta) {

    // MOTIVI DI USCITA DAL CICLO FOR

    //Per la ricerca lungo la retta
    if (ricerca_retta)
    {
        //cout << endl << "esco per la retta" << endl;
        return true;
    }

    //Solo se sono dopo il quarto livello (oppure dopo il secondo se non sono più al primo ciclo), allora posso terminare il 'cicle_programme' se il chi_quadro è migliorato più del per mille (in ogni caso questo non basta per terminare l'intero programma)
    if (((cicle > 1) || (livello > 4)) && ((chi_quadro.size() > 1) && (livello > 2)))
    {
        double check = (chi_quadro[chi_quadro.size() - 2] - chi_quadro[chi_quadro.size() - 1]) - chi_quadro[chi_quadro.size() - 1] * 0.001;
        if (check < 0) {
            //cout << endl << "esco chi quadro < 0.001" << endl;
            return true;
        }
    }

    //Se è stato migliorato solo una volta e sono già al 6° tentativo tanto vale uscire
    if (chi_quadro.size() == 1 && k > 5) {
        //cout << endl << "esco al sesto tentativo" << endl;
        return true;
    }

    //Se non è migliorato nemmeno una volta il chi quadro e sono al 7° tentativo tanto vale terminare tutto il programma
    if (chi_quadro.size() == 0 && k > 6) {
        return true;     //torna nella main principale uscendo anche dal while
    }   // con la stessa condizione esco anche con 'end()' [*]

    return false;
}

int covering::get_chi_quadro_size() {
    return chi_quadro.size();
}

bool covering::end() {

    //riaggiorno parametri, libero 'chi_quadro' e 'livelli', e passo al ciclo di programma successivo
    par = par_best;
    livelli.clear();            //livelli viene ripristinato
    livelli.assign(100, 0);     //
    chi_quadro.clear();         //chi_quadro viene ripristinato

    //Termino l'intero programma se il chi quadro subisce effettivamente un miglioramento minore del per mille (0.001)
    if (cicle_programms > 1 && (chi_quadro_ciclo_prec != chi_quadro_min && chi_quadro_ciclo_prec - chi_quadro_min < 0.001))
        return true;

    // [*] Se non è migliorato nemmeno una volta il chi quadro e ero al 7° tentativo tanto vale terminare tutto il programma
    if (chi_quadro.size() == 0 && k > 6) {
        return true;
    }

    return false;

}

