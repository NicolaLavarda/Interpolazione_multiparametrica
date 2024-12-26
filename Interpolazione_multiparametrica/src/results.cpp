#include "chi_square.h"
#include "matrix.h"
#include "interpolating_function.h"
#include "results.h"
#include "covering.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;

extern vector<double> x;        //Dati iniziali
extern double chi_quadro_min;   //Chi quadro minimo


// Chiamare nel programma 'Results(par_best, approx);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


Results_base::Results_base(vector<double> par_base, bool approx_bool, std::ostream& output) :
    par(par_base), chi_min(f_chi_quadro(par)), approx_bool(approx_bool) , out(output)
{
    double sensibility = 0.01;
    discesa_gradiente(par_base, sensibility);

    par = par_base;
    chi_min = f_chi_quadro(par_base);

}

void Results_base::general_result(double valore, double errore, string nome, bool arrotondamento) {
    double err_percentuale = fabs(errore / valore) * 100;
    if (!arrotondamento)
    {
        double min = (errore < fabs(valore)) ? errore : fabs(valore);
        double format = (min < 1e-6) ? ios::scientific : ios::fixed;
        out << nome << " = ( " << (format == ios::scientific ? scientific : fixed) << setprecision((format == ios::scientific) ? 6 : 10) << valore << " +- " << setprecision((format == ios::scientific) ? 6 : 10) << errore << " )";
        out << setw(8) << " " << "e_r = " << fixed << setprecision(10) << err_percentuale << "\%" << endl;
    }
    else
    {
        // Includere:   #include <iomanip>
        int digits_val = 0 - floor(log10(errore));
        if (errore >= 1) { digits_val = 0; }
        out << nome << " = ( " << fixed << setprecision(digits_val) << valore << " +- " << setprecision(digits_val) << errore << " )";

        int digits_e_r = 0 - floor(log10(err_percentuale));
        if (err_percentuale >= 1) digits_e_r = 0;
        out << setw(8) << " " << "e_r = " << fixed << setprecision(digits_e_r + 1) << err_percentuale << "\%" << endl;
    }
}


Result1::Result1()
    : Results_base(par, approx_bool, out) {
    //Stampo a schermo i risultati con errori dal chi+1
    out << endl;
    out << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    double chi_piu_uno_val = 1; if (false) chi_piu_uno_val = chi_quadro_piu_uno(par.size());
    out << "Parametri con errori al chi_quadro+1:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        string nome = "par" + to_string(i);
        double range = 0.20;   //cerco l'errore al chi+1 entro +-20% del valore del parametro (in caso fosse più grande la funzione allarga la ricerca)
        general_result(par[i], chi_piu_uno_par_n(par, i, chi_piu_uno_val, range), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    out << "----------------------------------------------------------------" << endl;
}


Result2::Result2()
    : Results_base(par, approx_bool, out) {
    //Calcolo Hessiana
    vector<vector<double>> H = hessian(par);
    //out << "Matrice Hessiana: " << endl;
    //stampaMatrice(H, out);

    //Calcolo matrice covarianza (inversa dell'hessiana)
    vector<vector<double>> invH = inversa(H);
    out << "Matrice di Covarianza: " << endl;
    stampaMatrice(invH, out);

    out << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    out << "Parametri con errori da matrice di covarianza:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        string nome = "par" + to_string(i);
        general_result(par[i], sqrt(invH[i][i]), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    out << "----------------------------------------------------------------" << endl;
}


Results::Results(std::vector<double>& par_derived, bool approx_bool, std::ostream& output) :
    Results_base(par_derived, approx_bool, output), Result1(), Result2() {
    //chiama in ordine 'Results_base', 'Result1' e 'Result2' ed infine qui dentro al costruttore di 'Results'
    // 
    // Chiamare nel programma 'Results(par_best, approx);'
    // in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


    //Stampo anche le informazioni sul chi_quadro
    double GDL = x.size() - par.size();
    double p_value_val = p_value(chi_min, GDL);
    out << "chi_quadro = " << setprecision(3) << chi_min << " / " << GDL << endl;
    out << "p_value = " << p_value_val << endl;
    if (p_value_val < 0.05)
        out << "Possibile sottostima degli errori" << endl;
    if (p_value_val > 0.95)
        out << "Possibile sovrastima degli errori" << endl;
    out << "----------------------------------------------------------------" << endl << endl;

    //Aggiornamento dei parametri all'esterno (se ho migliorato il chi_quadro e quindi i parametri)
    if (f_chi_quadro(par) < f_chi_quadro(par_derived))
    {
        par_derived = par;      // aggiorno al di fuori
        chi_quadro_min = chi_min;
    }

}

