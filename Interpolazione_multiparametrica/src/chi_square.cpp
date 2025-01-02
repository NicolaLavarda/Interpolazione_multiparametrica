#include "chi_square.h"

#include "covering.h"
#include "matrix.h"
#include "interpolator.h"

#include "results.h"
#include <iostream>
#include <vector>
#include <cmath>        // per std::nan() e altre cose
#include <limits>       // per std::numeric_limits

using namespace std;



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


double chi_piu_uno_par_n(vector<double> parametri, int numero_parametro, double chi_piu_uno_val, double range) {

    // Riferimento all'istanza Singleton
    Interpolator& i_generator = Interpolator::getInstance();

    //Definsco variabili iniziali in base a quale parametro (numero_parametro) voglio analizzare
    double par = parametri[numero_parametro];
    if (numero_parametro >= parametri.size()) {
        cout << endl << "ERROR chi+1 - not number parameters valid";
        cout << endl << endl; //exit(EXIT_FAILURE);
        return std::numeric_limits<double>::quiet_NaN();
    }

    //precisione
    double precisione_par = fabs(par * 0.0000001);
    double chi_min = i_generator.fChiQuadro(parametri);

    //Estremi per cui il sigma sarebbe altrimenti maggiore del range (di solito do valore 20%)
    range = (range == 0) ? 0.20 : range;
    double max_val_dx = par + fabs(par * range);              //non par*1.1 altrimenti errore con numeir negativi
    double max_val_sx = par - fabs(par * range);
    if (fabs(par) < 0.1)
    {
        max_val_dx = 1;       // se il parametro è vicino a 0 spesso l'errore può essere molto più grande del valore (ad esempio l'intercetta di una retta che viene 0.02 può avere benissimo un errore di 0.1)
        max_val_sx = -1;
    }

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
        double sum_chi = i_generator.fChiQuadro(parametri);
        dx_f_chi.push_back(parametri[numero_parametro]);
        dx_chi.push_back(sum_chi);
        //cout << parametri[numero_parametro] << "\t" << sum_chi << endl;
    }

    if (dx_chi[0] > dx_chi[1]) {
        posizione_min_chi_dx = 1;
        posizione_sec_min_chi_dx = 0;
    }


    int controllo_dx = 0;

    while (fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) > precisione_par)
    {

        parametri[numero_parametro] = dx_f_chi[posizione_min_chi_dx] + fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2;
        double sum_chi = i_generator.fChiQuadro(parametri);
        dx_f_chi.push_back(parametri[numero_parametro]);
        dx_chi.push_back(sum_chi);

        //cout << "->" << parametri[numero_parametro] << "\t" << sum_chi << endl;

        if (sum_chi > chi_min + chi_piu_uno_val)
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
            cout << endl << "ERROR chi+1 dx (par" << numero_parametro << ")";
            cout << endl << endl; //exit(EXIT_FAILURE);
            return std::numeric_limits<double>::quiet_NaN();
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
        double sum_chi = i_generator.fChiQuadro(parametri);
        sx_f_chi.push_back(parametri[numero_parametro]);
        sx_chi.push_back(sum_chi);
        //cout << parametri[numero_parametro] << endl;
    }


    int controllo_sx = 0;

    while (fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) > precisione_par)
    {

        parametri[numero_parametro] = sx_f_chi[posizione_sec_min_chi_sx] + fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2;
        //cout << parametri[numero_parametro] << endl;
        double sum_chi = i_generator.fChiQuadro(parametri);
        sx_f_chi.push_back(parametri[numero_parametro]);
        sx_chi.push_back(sum_chi);


        if (sum_chi > chi_min + chi_piu_uno_val)
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
            cout << endl << "| ERROR chi+1 sx (par" << numero_parametro << ")";
            cout << endl << endl; //exit(EXIT_FAILURE);
            return std::numeric_limits<double>::quiet_NaN();
            break;
        }

    }

    // Definisco il sigma finale, media dei due chi+1 trovati
    double a = dx_f_chi[posizione_min_chi_dx] + fabs(dx_f_chi[posizione_min_chi_dx] - dx_f_chi[posizione_sec_min_chi_dx]) / 2; //cout << endl << "-->" << a << endl;
    double b = sx_f_chi[posizione_min_chi_sx] + fabs(sx_f_chi[posizione_min_chi_sx] - sx_f_chi[posizione_sec_min_chi_sx]) / 2; //cout << endl << "-->" << b << endl;
    double sigma = (a - b) / 2;
    
    if (fabs(sigma / par) > 0.98 * range)       //se l'errore è al limite del range, è necessario ampliare il range
    {
        range *= 10;
        sigma = chi_piu_uno_par_n(parametri, numero_parametro, chi_piu_uno_val, range);     //continuo ricorsivamente la ricerca finché l'errore non sta dentor al range
    }

    return sigma;


}


