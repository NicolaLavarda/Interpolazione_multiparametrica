#include "interpolation.h"
#include "interpolating_function.h"
#include "matrix.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

extern vector<double> x, y, sigma_y;                             //Dati iniziali
extern double chi_quadro_min;                                    //Chi quadro minimo assoluto
extern vector<double> chi_quadro;                                //Vettore chi quadri minimi
extern vector<double> par_best;
extern vector<int> livelli;
extern vector<vector<double>> par_matrix;
extern int cicle_programms;                                                //numero di cicli di riduzione di 'spostamento'



//Ricerca automatica logaritmica
void ricerca_auto(vector<double>& par, vector<double>& par_auto, vector<double>& par_def, int n) {
    //mettere 'par_auto=par' e 'par_auto=par' prima della chiamata della funzione

    static double chi_min_auto = 1e10;
    static vector<double> sec_best_par = par;
    static double min = -1e10;
    static double max = 1e10;
    static double min_ordine = 0.1;
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
                while (n+k < par_size && par[n+k] != 0)
                    k++;
                //cout << "Sono a n = " << n << " e k = " << k << endl;
                k--;
            }
        }
        
        //if (salita)
        ricerca_auto(par, par_auto, par_def, n+1+k);

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


void algoritmo_bisezione(vector<double> par, vector<double>& par_def, const vector<double> passo, int n) {

    // Condizione di base: superato il numero di parametri da ottimizzare si esce dalla ricorsione (o sottoricorsione)
    if (n >= par.size())
        return;

    double range_min_par_n = par[n];
    double range_max_par_n = par[n];

    //Range con vettore passo
    if (passo.empty())
    {
        range_min_par_n -= fabs(par[n] * 0.2);        // Per la ricerca logaritmica automatica iniziale
        range_max_par_n += fabs(par[n] * 0.3);        //
    }
    else
    {
        range_min_par_n -= fabs(passo[n] / 2);        // Per lo svolgimento del programma principale
        range_max_par_n += fabs(passo[n] / 2);        //
    }

    //cout << range_min_par_n << " - " << par[n] << " - " << range_max_par_n << endl;

    double precisione_bisezione = 0.0001;
    if (par.size()<3)
        precisione_bisezione /= 100;    //Posso permettermi una precisione migliore con 1 o 2 parametri


    // PARAMETRO n
    double precisione_par_n = fabs(precisione_bisezione * par[n]);
    int controllo_par_n = 0;
    vector<double> chi_par_n, par_chi_n;
    double min_chi_par_n = 100000;                   //Importante valutar se mettere 'static' oppure no (meglio di no). Lascio che a volte sia "libera di sbagliare" e valutare combinazioni dei parametri che non necessariamente migliorino subito il chi quadro
    double sec_min_chi_par_n = 1000000;
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
            static int controllo_multiplo = 0;
            controllo_multiplo++;
            if (controllo_multiplo > 10) exit(EXIT_FAILURE);
            cout << endl << "ERRORE bisezione iniziale (par " << n << ")";
            cout << endl; //exit(EXIT_FAILURE);
            return;
            break;
        }
    }
}


// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void ricoprimento(vector<double>& par, vector<double>& par_def, const vector<double> passo, int livello, int dimensione, bool is_on_surface) {
    vector<double> par_prov = par_matrix[cicle_programms -1]; //non va bene 'static vector<double> par_prov = par;' perché altrimenti mi rimane anche per le riesecuzioni del programma successive
    
    if (dimensione == par.size()) {
        if (is_on_surface) {

            //for (int s = 0; s < par.size(); s++) { cout << "\t" << par[s]; } cout << endl;

            //Algoritmo di bisezione --------------------------------------------------------------------
            algoritmo_bisezione(par, par_def, passo, 0);

            double sum_chi_p = f_chi_quadro(par_def);

            

            if (sum_chi_p < chi_quadro_min) {
                if (1 == 2) {
                    string output_testo = "Parametri ottimizzati: \n";
                    for (int p = 0; p < par.size(); p++) {
                        output_testo += "par" + to_string(p) + ":\t" + to_string(par_def[p]) + "\n";
                    }
                    output_testo += "chi_quadro = " + to_string(sum_chi_p) + "\n";
                    cout << output_testo << endl;           //output valori
                }
                else
                    cout << "\t" << sum_chi_p;

                chi_quadro_min = sum_chi_p;
                par_best = par_def;
                //chi_quadro.push_back(sum_chi_p);          //Lo metto nella main in modo da aggiungerlo solo alla fine di un livello di ricoprimento (in uno stesso livello può essere che rimanga miniore del "per mille")
                livelli[livello]++;
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
        ricoprimento(par, par_def, passo, livello, dimensione + 1, is_on_surface || surface_condition);
    }
}


double chi_quadro_piu_uno(int num_parametri) {
    if (num_parametri == 1)
    {
        return 1;
    }
    else
    {
        double chi_critico = 0;
        double cum_prob = 0;
        double step_size = 0.0001;
        double x = 0;

        while (cum_prob < 0.682689492137086) {
            double pdf = exp(-0.5 * x) * pow(x, (num_parametri - 2) / 2.0) / (pow(2, num_parametri / 2.0) * tgamma(num_parametri / 2.0));
            chi_critico = x;
            cum_prob += pdf * step_size;
            x += step_size;
        }
        return chi_critico;
    }
}


double p_value(double chi_observato, int GDL) {
    if (GDL == 1) {
        return exp(-0.5 * chi_observato);  // Caso speciale per k=1
    }
    else {
        double cum_prob = 0;
        double step_size = 0.0001;
        double x = chi_observato;  // Inizia l'integrazione dal valore osservato
        double limite_superiore = max(100.0, chi_observato + 5.0 * GDL);
        double cum_prob_prec = 0;

        while (fabs(cum_prob - cum_prob_prec) > 0.0001 || x < limite_superiore) {  // Usa un limite massimo alto per approssimare l'infinito
            cum_prob_prec = cum_prob;
            double pdf = exp(-0.5 * x) * pow(x, (GDL - 2) / 2.0) /
                (pow(2, GDL / 2.0) * tgamma(GDL / 2.0));
            cum_prob += pdf * step_size;
            x += step_size;
        }

        return cum_prob;
    }
}


double chi_piu_uno_par_n(vector<double> parametri, int numero_parametro, int chi_piu_uno_val) {

    //Definsco variabili iniziali in base a quale parametro (numero_parametro) voglio analizzare
    double par = parametri[numero_parametro];
    if (numero_parametro >= parametri.size()) {
        cout << endl << "-----------------------------------";
        cout << endl << "| ERRORE ricerca errori parametri |";
        cout << endl << "-----------------------------------";
        cout << endl << endl; exit(EXIT_FAILURE);
    }

    //precisione
    double precisione_par = fabs(par * 0.0000001);

    //Estremi per cui il sigma sarebbe altrimenti maggiore del 30%
    double max_val_dx = par + fabs(par * 0.1);              //non par*1.1 altrimenti errore con numeir negativi
    double max_val_sx = par - fabs(par * 0.1);


    //Cerco il chi + 1 di destra --------------------------------

    vector<double> dx_f_chi, dx_chi;
    double min_chi_dx = par;
    double sec_min_chi_dx = max_val_dx;
    int posizione_min_chi_dx = 0;
    int posizione_sec_min_chi_dx = 1;


    vector<double> primi_chi_dx = { par, max_val_dx };
    for (int t = 0; t < primi_chi_dx.size(); t++)
    {
        parametri[numero_parametro] = primi_chi_dx[t];
        double sum_chi = f_chi_quadro(parametri);
        dx_f_chi.push_back(parametri[numero_parametro]);
        dx_chi.push_back(sum_chi);
        //cout << parametri[numero_parametro] << endl;
    }

    if (dx_chi[0] > dx_chi[1]) {
        posizione_min_chi_dx = 1;
        posizione_sec_min_chi_dx = 0;
    }


    int controllo_dx = 0;

    while (fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) > precisione_par)
    {

        parametri[numero_parametro] = dx_f_chi[posizione_min_chi_dx] + fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2;
        //cout << parametri[numero_parametro] << endl;
        double sum_chi = f_chi_quadro(parametri);
        dx_f_chi.push_back(parametri[numero_parametro]);
        dx_chi.push_back(sum_chi);

        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_dx = sum_chi;
            posizione_sec_min_chi_dx = dx_chi.size() - 1;
        }
        else {
            min_chi_dx = sum_chi;
            posizione_min_chi_dx = dx_chi.size() - 1;
        }

        controllo_dx++;
        if (controllo_dx > 100)
        {
            cout << endl << "---------------------------------";
            cout << endl << "| ERRORE chi+1 di destra (par" << numero_parametro << ") |";
            cout << endl << "---------------------------------";
            cout << endl << endl; exit(EXIT_FAILURE);
            break;
        }

    }
    //cout << endl << endl;

    //Cerco il chi + 1 di sinistra --------------------------------

    vector<double> sx_f_chi, sx_chi;
    double min_chi_sx = par;
    double sec_min_chi_sx = max_val_sx;
    int posizione_min_chi_sx = 0;
    int posizione_sec_min_chi_sx = 1;


    vector<double> primi_chi_sx = { par, max_val_sx };          // par è rimasto al valore originale (/ottimizzato)
    for (int t = 0; t < primi_chi_sx.size(); t++)
    {
        parametri[numero_parametro] = primi_chi_sx[t];
        double sum_chi = f_chi_quadro(parametri);
        sx_f_chi.push_back(parametri[numero_parametro]);
        sx_chi.push_back(sum_chi);
        //cout << parametri[numero_parametro] << endl;
    }


    int controllo_sx = 0;

    while (fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_par)
    {

        parametri[numero_parametro] = sx_f_chi[posizione_sec_min_chi_sx] + fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;
        //cout << parametri[numero_parametro] << endl;
        double sum_chi = f_chi_quadro(parametri);
        sx_f_chi.push_back(parametri[numero_parametro]);
        sx_chi.push_back(sum_chi);


        if (sum_chi > chi_quadro_min + chi_piu_uno_val)
        {
            sec_min_chi_sx = sum_chi;
            posizione_sec_min_chi_sx = sx_chi.size() - 1;
        }
        else {
            min_chi_sx = sum_chi;
            posizione_min_chi_sx = sx_chi.size() - 1;
        }

        controllo_sx++;
        if (controllo_sx > 100)
        {
            cout << endl << "-----------------------------------";
            cout << endl << "| ERRORE chi+1 di sinistra (par" << numero_parametro << ") |";
            cout << endl << "-----------------------------------";
            cout << endl << endl; exit(EXIT_FAILURE);
            break;
        }

    }

    // Definisco il sigma finale, media dei due chi+1 trovati
    double a = dx_f_chi[posizione_min_chi_dx] + fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2; //cout << endl << "-->" << a << endl;
    double b = sx_f_chi[posizione_min_chi_sx] + fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2; //cout << endl << "-->" << b << endl;
    return (a - b) / 2;
}




void linearFit(vector<double>& x, vector<double>& y, double& m, double& q) {
    // Controllo che i vettori abbiano la stessa dimensione e non siano vuoti
    if (x.size() != y.size() || x.empty()) {
        throw invalid_argument("I vettori x e y devono avere la stessa dimensione e non essere vuoti.");
    }

    int n = x.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;

    // Calcolo delle somme necessarie
    for (int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_xx += x[i] * x[i];
    }

    // Calcolo di m e q
    double denominator = n * sum_xx - sum_x * sum_x;
    if (denominator == 0) {
        throw runtime_error("Denominatore nullo. Interpolazione non definita.");
    }

    m = (n * sum_xy - sum_x * sum_y) / denominator;
    q = (sum_y - m * sum_x) / n;
}


void bisezione_lin(vector<double>& par, vector<double> m, vector<double> q) {
    
    //Prima parte a step abbastanza piccoli e tutti regolari
    double range_min_par_n = par[0] / 2;
    double range_max_par_n = par[0] * 10;

    if (par[0] < 0)
    {
        range_min_par_n = par[0] * 10;
        range_max_par_n = par[0] / 2;
    }

    //cout << range_min_par_n << " < " << par[0] << " < " << range_max_par_n << endl;
    
    //Blocco ricerca bisezione lungo la retta (ad esempio se par[0]=10 allora va da 5 a 100 in bisezione)
    if (1 == 1) {

        double precisione_bisezione = 0.00001;

        int controllo_par_n = 0;
        vector<double> chi_par_n, par_chi_n;
        double min_chi_par_n = 100000;                   //Importante valutar se mettere 'static' oppure no (meglio di no). Lascio che a volte sia "libera di sbagliare" e valutare combinazioni dei parametri che non necessariamente migliorino subito il chi quadro
        double sec_min_chi_par_n = 1000000;
        double posizione_min_chi_par_n = 0;
        double posizione_sec_min_chi_par_n = 1;

        for (double p = range_min_par_n; p < range_max_par_n + fabs(range_max_par_n * 0.1); p += (range_max_par_n - range_min_par_n))   //primi due elementi agli estremi --> uso 'p < range_max_par_n + fabs(range_max_par_n * 0.1)' anzichè 'p <= range_max_par_n' per evitare di controllare valori tra double molto vicini
        {
            par[0] = p;
            // Definisco gli altri parametri in funzione del primo ( 'par[0]' ) utilizzando gli 'm' e 'q' secondo par[n] = m_n * par[0] + q_n;
            for (int i = 1; i < par.size(); i++)
            {
                par[i] = m[i - 1] * par[0] + q[i - 1];
            }
            double sum_chi = f_chi_quadro(par);
            chi_par_n.push_back(sum_chi);
            par_chi_n.push_back(p);
        }
        if (chi_par_n[0] > chi_par_n[1]) {
            posizione_min_chi_par_n = 1;
            posizione_sec_min_chi_par_n = 0;
        }

        while (fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) > precisione_bisezione)
        {

            if (par_chi_n[posizione_min_chi_par_n] > par_chi_n[posizione_sec_min_chi_par_n])
            {
                par[0] = par_chi_n[posizione_sec_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
            }
            else
            {
                par[0] = par_chi_n[posizione_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
            }

            // Definisco gli altri parametri in funzione del primo ( 'par[0]' ) utilizzando gli 'm' e 'q' secondo par[n] = m_n * par[0] + q_n;
            for (int i = 1; i < par.size(); i++)
            {
                par[i] = m[i - 1] * par[0] + q[i - 1];
            }
            double sum_chi = f_chi_quadro(par);
            par_chi_n.push_back(par[0]);
            chi_par_n.push_back(sum_chi);

            /*
            cout << endl;
            for (int k = 0; k < par.size(); k++)
                cout << par[k] << "\t";
            cout << " -> ";
            cout << sum_chi << endl;
            */


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
                static int controllo_multiplo_lin = 0;
                controllo_multiplo_lin++;
                if (controllo_multiplo_lin > 10) exit(EXIT_FAILURE);
                cout << endl << "ERRORE bisezione interpolazione lineare";
                cout << endl; //exit(EXIT_FAILURE);
                return;
                break;
            }
        }

    }


    /*
    cout << endl;
    for (int k = 0; k < par.size(); k++)
        cout << par[k] << "\t";
    cout << " -> ";
    cout << f_chi_quadro(par) << endl;
    */


    // par[0] ora è stato migliorato con "Blocco ricerca bisezione lungo la retta"

    range_min_par_n = par[0] * 0.7;
    range_max_par_n = par[0] * 1.3;

    if (par[0] < 0)
    {
        range_min_par_n = par[0] * 1.3;
        range_max_par_n = par[0] * 0.7;
    }


    double min_chi_lin_passi = 1e10;

    vector<double> par_lin_finali;
    double passo_i = (range_max_par_n - range_min_par_n) / 50;
    for (double p = range_min_par_n; p < range_max_par_n; p += passo_i)
    {
        vector<double> par_lin_passi = par;
        vector<double> passi_lin;
        par_lin_passi[0] = p;
        passi_lin.push_back(passo_i*2);
        for (int i = 1; i < par.size(); i++)
        {
            par_lin_passi[i] = m[i - 1] * par_lin_passi[0] + q[i - 1];
            passi_lin.push_back(fabs(par_lin_passi[i]));                    //Il passo (che poi diventa la metà in bisezione) è pari al parametro così va circa par*0.5 < par < par*1.5
        }
        /*
        cout << "-------" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par_lin_passi[k] << "\t";
        cout << " -> ";
        cout << f_chi_quadro(par_lin_passi) << endl;
        */

        /*
        vector<double> quadrante1 = par_lin_passi;
        vector<double> quadrante2 = par_lin_passi;
        vector<double> quadrante3 = par_lin_passi;
        vector<double> quadrante4 = par_lin_passi;

        vector<double> passi1;
        vector<double> passi2;
        vector<double> passi3;
        vector<double> passi4;

        for (int i = 1; i < par.size(); i++)        //Perpendicolarmente alla retta di direzione di minimizzazione del chi quardo, cerco su una superficie n-1 dim con bisezione di migliorare gli altri n-1 parametri (divido lo spazio in 4 quadranti con un vertice al centro in comune sul punto in questione della retta)
        {
            quadrante1[i] *= 1.5;
            quadrante2[i] *= 1.5;
            quadrante3[i] *= 1.5;
            quadrante4[i] *= 1.5;
        }
        */

        algoritmo_bisezione(par_lin_passi, par_lin_passi, passi_lin, 0);

        /*
        cout << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par_lin_passi[k] << "\t";
        cout << " -> ";
        cout << f_chi_quadro(par_lin_passi) << endl;
        cout << "-------" << endl;
        */

        double sum_chi_lin = f_chi_quadro(par_lin_passi);
        if (sum_chi_lin < min_chi_lin_passi)
        {
            min_chi_lin_passi = sum_chi_lin;
            par_lin_finali = par_lin_passi;
        }
    }
    par = par_lin_finali;
    
    /*
    cout << endl;
    for (int k = 0; k < par.size(); k++)
        cout << par[k] << "\t";
    cout << " -> ";
    cout << f_chi_quadro(par) << endl;
    */

}









void risultato(double valore, double errore, string nome, bool arrotondamento) {
    if (!arrotondamento)
    {
        cout << nome << " = ( " << fixed << setprecision(10) << valore << " +- " << setprecision(10) << errore << " )" << endl;
    }
    else
    {
        // Includere:   #include <iomanip>
        int digits_val = 0 - floor(log10(errore));
        if (errore >= 1) { digits_val = 0; }
        cout << nome << " = ( " << fixed << setprecision(digits_val) << valore << " +- " << setprecision(digits_val) << errore << " )";
    }

    double err_percentuale = fabs(errore / valore) * 100;    
    if (!arrotondamento)
    {
        cout << setw(12) << "\t" << "e_r = " << fixed << setprecision(10) << err_percentuale << "\%" << endl;
    }
    else
    {
        int digits_e_r = 0 - floor(log10(err_percentuale));
        if (err_percentuale >= 1) digits_e_r = 0;
        cout << setw(12) << "\t" << "e_r = " << fixed << setprecision(digits_e_r + 1) << err_percentuale << "\%" << endl;
    }
}

//Stampo a schermo i risultati con errori dal chi+1
void risultato1(vector<double> par_best, bool approx_bool) {
    cout << endl << endl;
    cout << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    double chi_piu_uno_val = 1; if (false) chi_piu_uno_val = chi_quadro_piu_uno(par_best.size());
    cout << "Parametri con errori al chi_quadro+1:" << endl;
    for (int i = 0; i < par_best.size(); i++)
    {
        string nome = "par" + to_string(i);
        risultato(par_best[i], chi_piu_uno_par_n(par_best, i, chi_piu_uno_val), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    cout << "----------------------------------------------------------------" << endl;
}

//Stampo a schermo i risultati con errori da matrice di covarianza
void risultato2(vector<double> par_best, bool approx_bool) {
    //Calcolo Hessiana
    vector<vector<double>> H = hessian(par_best);
    //cout << "Matrice Hessiana: " << endl;
    //stampaMatrice(H);


    //Calcolo matrice covarianza (inversa dell'hessiana)
    vector<vector<double>> invH = inversa(H);
    cout << "Matirce di Covarianza: " << endl;
    stampaMatrice(invH);


    cout << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    cout << "Parametri con errori da matrice di covarianza:" << endl;
    for (int i = 0; i < par_best.size(); i++)
    {
        string nome = "par" + to_string(i);
        risultato(par_best[i], sqrt(invH[i][i]), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    double GDL = x.size() - par_best.size();
    double p_value_val = p_value(chi_quadro_min, GDL);
    cout << "chi_quadro = " << setprecision(3) << f_chi_quadro(par_best) << " / " << GDL << endl;
    cout << "p_value = " << p_value_val << endl;
    if (p_value_val < 0.05)
        cout << "Possibile sottostima degli errori" << endl;
    if (p_value_val > 0.95)
        cout << "Possibile sovrastima degli errori" << endl;

    cout << "----------------------------------------------------------------" << endl << endl;
}

