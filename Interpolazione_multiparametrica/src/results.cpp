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

using namespace std;

extern vector<double> x;        //Dati iniziali


// Chiamare nel programma 'Results(par_best, approx);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


Results_base::Results_base(vector<double> par_base, bool approx_bool) :
    par(par_base), chi_min(f_chi_quadro(par)), approx_bool(approx_bool)
{
    //cout << endl << "Risultati: ";

    double sensibility = 0.01;
    discesa_gradiente(par_base, sensibility);
    
    par = par_base;
    chi_min = f_chi_quadro(par_base);

}

void Results_base::general_result(double valore, double errore, string nome, bool arrotondamento) {
    double err_percentuale = fabs(errore / valore) * 100;
    if (!arrotondamento)
    {
        cout << nome << " = ( " << fixed << setprecision(10) << valore << " +- " << setprecision(10) << errore << " )";
        cout << setw(8) << " " << "e_r = " << fixed << setprecision(10) << err_percentuale << "\%" << endl;
    }
    else
    {
        // Includere:   #include <iomanip>
        int digits_val = 0 - floor(log10(errore));
        if (errore >= 1) { digits_val = 0; }
        cout << nome << " = ( " << fixed << setprecision(digits_val) << valore << " +- " << setprecision(digits_val) << errore << " )";

        int digits_e_r = 0 - floor(log10(err_percentuale));
        if (err_percentuale >= 1) digits_e_r = 0;
        cout << setw(8) << " " << "e_r = " << fixed << setprecision(digits_e_r + 1) << err_percentuale << "\%" << endl;
    }
}


Result1::Result1()
    : Results_base(par, approx_bool) {
    //Stampo a schermo i risultati con errori dal chi+1
    cout << endl;
    cout << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    double chi_piu_uno_val = 1; if (false) chi_piu_uno_val = chi_quadro_piu_uno(par.size());
    cout << "Parametri con errori al chi_quadro+1:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        string nome = "par" + to_string(i);
        general_result(par[i], chi_piu_uno_par_n(par, i, chi_piu_uno_val), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    cout << "----------------------------------------------------------------" << endl;
}


Result2::Result2()
    : Results_base(par, approx_bool) {
    //Calcolo Hessiana
    vector<vector<double>> H = hessian(par);
    //cout << "Matrice Hessiana: " << endl;
    //stampaMatrice(H);

    //Calcolo matrice covarianza (inversa dell'hessiana)
    vector<vector<double>> invH = inversa(H);
    cout << "Matrice di Covarianza: " << endl;
    stampaMatrice(invH);

    cout << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    cout << "Parametri con errori da matrice di covarianza:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        string nome = "par" + to_string(i);
        general_result(par[i], sqrt(invH[i][i]), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    cout << "----------------------------------------------------------------" << endl;
}


Results::Results(std::vector<double>& par_derived, bool approx_bool) :
    Results_base(par_derived, approx_bool), Result1(), Result2() {
    //chiama in ordine 'Results_base', 'Result1' e 'Result2' ed infine qui dentro al costruttore di 'Results'
    // 
    // Chiamare nel programma 'Results(par_best, approx);'
    // in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati



    //Stampo anche le informazioni sul chi_quadro
    double GDL = x.size() - par.size();
    double p_value_val = p_value(chi_min, GDL);
    cout << "chi_quadro = " << setprecision(3) << chi_min << " / " << GDL << endl;
    cout << "p_value = " << p_value_val << endl;
    if (p_value_val < 0.05)
        cout << "Possibile sottostima degli errori" << endl;
    if (p_value_val > 0.95)
        cout << "Possibile sovrastima degli errori" << endl;
    cout << "----------------------------------------------------------------" << endl << endl;


    //Aggiornamento dei parametri all'esterno (se ho migliorato il chi_quadro e quindi i parametri)
    if (f_chi_quadro(par) < f_chi_quadro(par_derived))
        par_derived = par;      // aggiorno al di fuori
}

