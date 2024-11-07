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

using namespace std;

extern vector<double> x, y, sigma_y;                             //Dati iniziali
extern double chi_quadro_min;                                    //Chi quadro minimo assoluto
extern vector<double> chi_quadro;                                //Vettore chi quadri minimi
extern vector<double> par_best;
extern vector<int> livelli;
extern vector<vector<double>> par_matrix;
extern int cicle_programms;                                                //numero di cicli di riduzione di 'spostamento'


//Ricerca automatica logaritmica
void ricerca_auto(vector<double>& par, vector<double>& par_auto, int n) {
    if (n >= par.size()) return;
    static double chi_min_auto = 1e10;
    double min = -1e10;
    double max = 1e10;
    double a = min;
    while (a <= max)
    {
        par_auto[n] = a;
        ricerca_auto(par, par_auto, n + 1);

        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par_auto, i)) / sigma_y[i], 2);
        }

        /*
        for (int k = 0; k < par.size(); k++)
        {
            cout << par_auto[k] << "\t";
        }
        cout << endl;
        cout << sum_chi << endl;
        */

        if (sum_chi < chi_min_auto)
        {
            par = par_auto;
            chi_min_auto = sum_chi;
        }

        // Devo andare da -1e9 a 1e9 togliendo [-1e-2;+1e-2]
        (a < 0) ? a /= 10 : a *= 10;
        if (a >= -0.001 - 1e-9 && a <= -0.001 + 1e-9)
            a = 0.01;
    }
}



void algoritmo_bisezione(vector<double> par, vector<double>& par_def, const vector<double> passo, int n) {

    // Condizione di base: superato il numero di parametri da ottimizzare si esce dalla ricorsione (o sottoricorsione)
    if (n >= par.size())
        return;

    //Range con vettore passo
    double range_min_par_n = par[n] - passo[n] / 2;
    double range_max_par_n = par[n] + passo[n] / 2;
    

    double precisione_bisezione = 0.00001;


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
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
        }
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
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
        }
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
            cout << endl << "-------------------------------------";
            cout << endl << "| ERRORE bisezione iniziale (par n) |";
            cout << endl << "-------------------------------------";
            cout << endl << endl; //exit(EXIT_FAILURE);
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

            double sum_chi_p = 0.0;
            for (int i = 0; i < x.size(); i++)
            {
                sum_chi_p += pow((y[i] - funzione_interpolante(x, par_def, i)) / sigma_y[i], 2);
            }

            

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
                chi_quadro.push_back(sum_chi_p);
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
        double sum_chi = 0;
        parametri[numero_parametro] = primi_chi_dx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, parametri, i)) / sigma_y[i], 2);
        }
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
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, parametri, i)) / sigma_y[i], 2);
        }
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
        double sum_chi = 0;
        parametri[numero_parametro] = primi_chi_sx[t];
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, parametri, i)) / sigma_y[i], 2);
        }
        sx_f_chi.push_back(parametri[numero_parametro]);
        sx_chi.push_back(sum_chi);
        //cout << parametri[numero_parametro] << endl;
    }


    int controllo_sx = 0;

    while (fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_par)
    {

        parametri[numero_parametro] = sx_f_chi[posizione_sec_min_chi_sx] + fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;
        //cout << parametri[numero_parametro] << endl;
        double sum_chi = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chi += pow((y[i] - funzione_interpolante(x, parametri, i)) / sigma_y[i], 2);
        }
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
    double err_percentuale = errore / valore * 100;
    int digits_e_r = 0 - floor(log10(err_percentuale));
    if (err_percentuale >= 1) { digits_e_r = 0; }
    cout << setw(12) << "\t" << "e_r = " << fixed << setprecision(digits_e_r + 1) << err_percentuale << "\%" << endl;
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
    cout << "Matrice Hessiana: " << endl;
    stampaMatrice(H);

    cout << endl;

    //Calcolo matrice covarianza (inversa dell'hessiana)
    vector<vector<double>> invH = inversa(H);
    cout << "Matirce di Covarianza: " << endl;
    stampaMatrice(invH);


    cout << endl;
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
    cout << "chi_quadro = " << setprecision(3) << chi_quadro_min << " / " << GDL << endl;
    cout << "p_value = " << p_value_val << endl;
    if (p_value_val < 0.05)
        cout << "Possibile sottostima degli errori" << endl;
    if (p_value_val > 0.95)
        cout << "Possibile sovrastima degli errori" << endl;

    cout << "----------------------------------------------------------------" << endl << endl;
}