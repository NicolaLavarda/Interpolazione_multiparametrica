#include "interpolating_function.h"
#include "bisection_algorithm.h"

#include <vector>
#include <cmath>

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>  // per std::max_element

using namespace std;


// Chiamare nel programma 'bisection_algorithm(par_best, step, chi_quadro_min);'
// in modo da creare un oggetto temporaneo che migliori i parametri con metodo di bisezione


bisection_algorithm::bisection_algorithm(vector<double>& par_best, const vector<double> step, double& chi_quadro_min) :
    par_base(par_best), step(step), chi_min(f_chi_quadro(par_best))
{

    //Preservo i parametri per controlalre alla fine se il chi quadro è stato effettivamente migliorato (ovvio, ma per doppio check)
    vector<double> par_prec = par_best;
    double chi_prec = f_chi_quadro(par_prec);
    
    try
    {
        bisection(par_best, par_best, step, 0);
    }
    catch (const std::runtime_error& e) {
        par_best = par_prec;            // in caso di errore ripristino 'par_best' e 'chi_quadro_min' e termino la funzione
        chi_quadro_min = chi_prec;      //
        std::cerr << "ERROR in bisection algorithm: " << e.what() << endl;
        return;
    }

    //Verifico siano effettivamente migliori i parametri (chi_quadro minore)
    double chi_bisezione = f_chi_quadro(par_best);

    if (chi_bisezione < chi_prec)
    {
        //par sono già modificati e i migliori
        chi_quadro_min = chi_bisezione;
    }
    else
    {
        par_best = par_prec;                 //Riporto i parametri e il chi_quadro ai valori originali (prima della chiamata della funzione)
        chi_quadro_min = chi_prec;           //
    }
    //cout << par_best[0] << "\t" << par_best[1] << "\t" << chi_quadro_min << endl;

    return;
}

void bisection_algorithm::bisection(vector<double> par, vector<double>& par_def, const vector<double> passo, int n) {

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

        if (par[n] > 0)   // 3.16 sta logaritmicamente a metà tra 1 e 10
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
        bisection(par, par_def, passo, n + 1);

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

