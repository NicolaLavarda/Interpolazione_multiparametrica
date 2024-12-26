#include "covering.h"
#include "interpolating_function.h"
#include "matrix.h"
//#include "results.h"
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



void algoritmo_bisezione(vector<double> par, vector<double>& par_def, const vector<double> passo, int n) {

    // Condizione di base: superato il numero di parametri da ottimizzare si esce dalla ricorsione (o sottoricorsione)
    if (n >= par.size())
        return;

    double range_min_par_n = par[n];
    double range_max_par_n = par[n];

    //Range con vettore passo
    if (passo.empty())
    {
        //range_min_par_n -= fabs(par[n] * 0.2);        // Per la ricerca logaritmica automatica iniziale
        //range_max_par_n += fabs(par[n] * 0.3);        //

        if (par[n]>0)   // 3.16 sta logaritmicamente a metà tra 1 e 10
        {
            range_min_par_n = par[n] * 0.316;       // ho cambiato rispetto a prima
            range_max_par_n = par[n] * 3.16;        //
        }
        else
        {
            range_min_par_n = par[n] * 3.16;        // ho cambiato rispetto a prima
            range_max_par_n = par[n] * 0.316;       //
        }
    }
    else
    {
        range_min_par_n -= fabs(passo[n] / 2);        // Per lo svolgimento del programma principale
        range_max_par_n += fabs(passo[n] / 2);        //
    }

    //cout << range_min_par_n << " - " << par[n] << " - " << range_max_par_n << endl;

    double precisione_bisezione = 0.0001;       //es. 0.0001 vuol dire che miglioro ogni parametro fino ad una parte su diecimila del valore del parametro
    if (par.size() < 3)
        precisione_bisezione /= 100;    //Posso permettermi una precisione migliore con 1 o 2 parametri


    // PARAMETRO n
    double precisione_par_n = fabs(precisione_bisezione * par[n]);
    int controllo_par_n = 0;
    vector<double> chi_par_n, par_chi_n;
    double min_chi_par_n = 1e20;                   //Importante valutar se mettere 'static' oppure no (meglio di no). Lascio che a volte sia "libera di sbagliare" e valutare combinazioni dei parametri che non necessariamente migliorino subito il chi quadro
    double sec_min_chi_par_n = 1e21;               //Questo non servirebbe nemmeno inizializzarlo
    double posizione_min_chi_par_n = 0;
    double posizione_sec_min_chi_par_n = 1;

    for (double p = range_min_par_n; p < range_max_par_n + fabs(range_max_par_n * 0.1); p += (range_max_par_n - range_min_par_n))   //primi due elementi agli estremi --> uso 'p < range_max_par_n + fabs(range_max_par_n * 0.1)' anzichè 'p <= range_max_par_n' per evitare di controllare valori tra double molto vicini
    {
        par[n] = p;
        double sum_chi = f_chi_quadro(par);
        chi_par_n.push_back(sum_chi);
        par_chi_n.push_back(p);
    }
    if (chi_par_n[0] > chi_par_n[1]) {
        posizione_min_chi_par_n = 1;
        posizione_sec_min_chi_par_n = 0;
    }

    while (fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) > precisione_par_n)
    {

        if (par_chi_n[posizione_min_chi_par_n] > par_chi_n[posizione_sec_min_chi_par_n])
        {
            par[n] = par_chi_n[posizione_sec_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
        }
        else
        {
            par[n] = par_chi_n[posizione_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
        }

        //Aggiornamento all'esterno degli ultimi parametri ottimizzzati
        par_def = par;

        // RICHIAMO ALGORITMO DI BISEZIONE RICORSIVO
        algoritmo_bisezione(par, par_def, passo, n + 1);

        // ritorno a parametro precedente
        double sum_chi = f_chi_quadro(par);
        par_chi_n.push_back(par[n]);
        chi_par_n.push_back(sum_chi);


        if (sum_chi < min_chi_par_n)
        {
            sec_min_chi_par_n = min_chi_par_n;
            min_chi_par_n = sum_chi;
            posizione_sec_min_chi_par_n = posizione_min_chi_par_n;
            posizione_min_chi_par_n = chi_par_n.size() - 1;
        }
        else {
            sec_min_chi_par_n = sum_chi;
            posizione_sec_min_chi_par_n = chi_par_n.size() - 1;
        }

        controllo_par_n++;
        if (controllo_par_n > 100)
        {
            throw runtime_error("bisection algorithm faild");
            //static int controllo_multiplo = 0;
            //controllo_multiplo++;
            //if (controllo_multiplo > 10) exit(EXIT_FAILURE);
            //cout << endl << "ERRORE bisezione iniziale (par " << n << ")";
            //cout << endl; //exit(EXIT_FAILURE);
            return;
            break;
        }
    }

}


// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void ricoprimento(vector<double>& par, vector<double>& par_def, const vector<double> passo, int livello, int dimensione, bool is_on_surface, bool output) {
    vector<double> par_prov = par_matrix[cicle_programms - 1]; //non va bene 'static vector<double> par_prov = par;' perché altrimenti mi rimane anche per le "riesecuzioni del programma successive". Importante che sia esattamente così: se passassi in qualsiasi modo il vettore alla funzione fidati che non funziona perché poi cambia, fidati.

    if (dimensione == par.size()) {
        if (is_on_surface) {

            //for (int s = 0; s < par.size(); s++) { cout << "\t" << par[s]; } cout << endl;

            //Algoritmo di bisezione --------------------------------------------------------------------
            vector<double> par_provv_def = par_def;
            try
            {
                algoritmo_bisezione(par, par_def, passo, 0);
            }
            catch (const std::runtime_error& e) {
                par_def = par_provv_def;            // la funzione di bisezione poteva modificare solo 'par_def' e in caso di errore lo ripristino
                std::cerr << "ERROR in bisection algorithm: " << e.what() << endl;
            }

            double sum_chi_p = f_chi_quadro(par_def);



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
                {
                    cout << "\t" << sum_chi_p;
                }

                chi_quadro_min = sum_chi_p;
                par_best = par_def;
                chi_quadro.push_back(sum_chi_p);          //Lo metto nella main in modo da aggiungerlo solo alla fine di un livello di ricoprimento (in uno stesso livello può essere che rimanga miniore del "per mille")
                livelli[livello]++;                         //Metto un valore diverso da 0 se in un livello riesco a migliorare il chi quadro
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


// gradiente della funzione chi qaudro (h_gen% del valore)
std::vector<double> grad_f_chi_quadro(const std::vector<double> par) {
    std::vector<double> grad(par.size());
    double h = 1e-3;
    //cout << endl << "calcolo gradiente:" << endl;
    // Calcoliamo la derivata parziale per ogni parametro
    for (size_t i = 0; i < par.size(); ++i) {
        // Crea una copia del vettore dei parametri
        std::vector<double> par_plus = par;
        std::vector<double> par_minus = par;

        // Modifica il parametro i-esimo in positivo e negativo
        par_plus[i] += h;
        par_minus[i] -= h;

        // Approssima la derivata con la differenza centrale
        grad[i] = (f_chi_quadro(par_plus) - f_chi_quadro(par_minus)) / (2 * h);
        //cout <<par_plus[i] << "\t" << par_minus[i] << "\t" << f_chi_quadro(par_plus) << "\t" << f_chi_quadro(par_minus) << "\t" << grad[i] << endl;
    }

    return grad;
}

// Metodo di discesa del gradiente
void discesa_gradiente(std::vector<double>& par, double sensibility) {

    double tasso_apprendimento = sensibility;        //es. sensibility=0.1 -> all'inizio a passi del 10% dei valori
    double tollerance = 0.001;
    int max_iter = 300;

    double grad0 = f_chi_quadro(par);
    double grad1 = grad0;
    bool migliorato = false;

    for (int iter = 0; iter < max_iter; ++iter) {

        std::vector<double> grad = grad_f_chi_quadro(par);  // Calcola il gradiente

        std::vector<double> check_max;
        for (int i = 0; i < par.size(); i++)
        {
            check_max.push_back(fabs(grad[i] / par[i]));
            //cout << " -> " << fabs(grad[i] / par[i]);
        }
        //cout << endl;
        auto max_iter = std::max_element(check_max.begin(), check_max.end());
        int max_index = std::distance(check_max.begin(), max_iter);
        //cout << " -> " << max_index << endl;

        //cout << "gradiente: \t";
        // Aggiorna i parametri con il gradiente
        vector<double> par_p = par;
        for (size_t i = 0; i < par.size(); ++i) {
            par_p[i] -= grad[i] * par_p[i] * tasso_apprendimento / fabs(grad[max_index]);     // '/ grad[max_index]' normalizza in modo che non sia troppo grande ad esempio se par[0]=12 e grad[0]=180 in questo modo viene 180*(12*0.1)/180=1.2 mentre un par[1]=80 con grad[1]=5 diminuirà di 5*(80*0.1)/180=0.2
            //cout << grad[i] * par_p[i] * tasso_apprendimento / fabs(grad[max_index]) << " - " << par_p[i] << "\t";
        }

        double grad_p = f_chi_quadro(par_p);
        if (grad_p < grad0)
        {
            grad1 = grad0;
            grad0 = grad_p;
            migliorato = true;
            par = par_p;
        }
        else
        {
            tasso_apprendimento /= 2;
        }
        //cout << migliorato << endl;
        //std::cout << "Iterazione " << iter + 1 << ": f_chi_quadro(par) = " << f_chi_quadro(par) << std::endl;

        // Condizione di convergenza
        if (migliorato && (std::fabs(grad0 - grad1) < tollerance)) {
            //std::cout << "Convergenza raggiunta dopo " << iter + 1 << " iterazioni." << std::endl;
            break;
        }
    }
}