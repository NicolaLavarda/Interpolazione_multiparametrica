#include "results.h"

#include "interpolator.h"
#include "chi_square.h"
#include "matrix.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>

extern int dim_x;


// Chiamare nel programma 'Results(par_best, approx);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


Results_base::Results_base(std::vector<double> par_base, bool approx_bool, std::ostream& output) :
    par(par_base), approx_bool(approx_bool) , out(output)
{
    chi_min = i_generator.fChiQuadro(par);
}

void Results_base::general_result(double valore, double errore, std::string nome, bool arrotondamento) {
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
    out << "Parameters with errors at chi-squared+1:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        std::string nome(1, name_par[i]);
        double range = 0.20;   //cerco l'errore al chi+1 entro +-20% del valore del parametro (in caso fosse più grande la funzione allarga la ricerca)
        general_result(par[i], chi_piu_uno_par_n(par, i, chi_piu_uno_val, range), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    out << "----------------------------------------------------------------" << endl;
}


Result2::Result2()
    : Results_base(par, approx_bool, out) {
    //Calcolo Hessiana
    std::vector<std::vector<double>> H = hessian(par);
    //out << "Matrice Hessiana: " << endl;
    //stampaMatrice(H, out);

    //Calcolo matrice covarianza (inversa dell'hessiana)
    std::vector<std::vector<double>> invH = inversa(H);
    out << "Covariance Matrix: " << endl;
    stampaMatrice(invH, out);

    out << "----------------------------------------------------------------" << endl;
    //Calcolo errori e stampo a schermo i risultati
    out << "Parameters with errors from the covariance matrix:" << endl;
    for (int i = 0; i < par.size(); i++)
    {
        std::string nome(1, name_par[i]);
        general_result(par[i], sqrt(invH[i][i]), nome, approx_bool);      //true per arrotondare con il giusto numero di cifre significative
    }
    out << "----------------------------------------------------------------" << endl;
}


Results::Results(std::vector<double> par_derived, bool approx_bool, int dim_x, std::ostream& output) :
    Results_base(par_derived, approx_bool, output), Result1(), Result2() {
    //chiama in ordine 'Results_base', 'Result1' e 'Result2' ed infine qui dentro al costruttore di 'Results'
    // 
    // Chiamare nel programma 'Results(par_best, approx);'
    // in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


    //Stampo anche le informazioni sul chi_quadro
    double GDL = dim_x - par.size();
    double p_value_val = p_value(chi_min, GDL);
    out << "chi-squared = " << setprecision(3) << chi_min << " / " << GDL << endl;
    out << "p_value = " << p_value_val << endl;
    if (p_value_val < 0.05)
        out << "Possible underestimation of errors" << endl;
    if (p_value_val > 0.95)
        out << "Possible overestimation of errors" << endl;
    out << "----------------------------------------------------------------" << endl << endl;

}

