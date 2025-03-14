#include "interpolator.h"
#include "gradient_descent_algorithm.h"

#include <vector>
#include <cmath>

#include <iostream>
#include <cstdlib>
#include <algorithm>  // per std::max_element



// Chiamare nel programma 'gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);'
// in modo da creare un oggetto temporaneo che migliori i parametri con metodo di bisezione


gradient_descent_algorithm::gradient_descent_algorithm(std::vector<double>& par, double& chi_quadro_min, double sensibility)
{
    //Preservo i parametri per controlalre alla fine se il chi quadro è stato effettivamente migliorato (ovvio, ma per doppio check)
    std::vector<double> par_prec = par;

    double tasso_apprendimento = sensibility;        //es. sensibility=0.1 -> all'inizio a passi del 10% dei valori
    double tollerance = 0.001;

    double grad0 = i_generator->fChiQuadro(par);
    double grad1 = grad0;
    bool migliorato = false;

    int par_size = par.size();

    for (int iter = 0; iter < 300; ++iter) {

        std::vector<double> grad = grad_f_chi_quadro(par);  // Calcola il gradiente

        std::vector<double> check_max;
        for (int i = 0; i < par_size; i++)
        {
            check_max.push_back(fabs(grad[i] / par[i]));
            //cout << " -> " << fabs(grad[i] / par[i]);
        }
        //cout << endl;
        auto max_iter = std::max_element(check_max.begin(), check_max.end());
        int max_index = std::distance(check_max.begin(), max_iter);
        //cout << " -> " << max_index << endl;

        //std::cout << "gradiente: \t";
        // Aggiorna i parametri con il gradiente
        std::vector<double> par_p = par;
        for (int i = 0; i < par_size; ++i) {
            par_p[i] -= grad[i] * par_p[i] * tasso_apprendimento / fabs(grad[max_index]);     // '/ grad[max_index]' normalizza in modo che non sia troppo grande ad esempio se par[0]=12 e grad[0]=180 in questo modo viene 180*(12*0.1)/180=1.2 mentre un par[1]=80 con grad[1]=5 diminuirà di 5*(80*0.1)/180=0.2
            //std::cout << grad[i] * par_p[i] * tasso_apprendimento / fabs(grad[max_index]) << " - " << par_p[i] << "\t";
        }

        double grad_p = i_generator->fChiQuadro(par_p); //std::cout << "|chi = " << grad_p << "| \t";
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
        //debug
        //std::cout << "Iterazione " << iter + 1 << ": i_generator->fChiQuadro(par) = " << i_generator->fChiQuadro(par); for (int i = 0; i < par_size; i++) std::cout << "\t" << par[i]; std::cout << std::endl;

        // Condizione di convergenza
        if (migliorato && (std::fabs(grad0 - grad1) < tollerance)) {
            //std::cout << "Convergenza raggiunta dopo " << iter + 1 << " iterazioni." << std::endl;
            break;
        }
    }

    /*
    for (size_t i = 0; i < par_size; i++)
        std::cout << par[i] << "\t";
    std::cout << std::endl;
    */

    //Verifico siano effettivamente migliori i parametri (chi_quadro minore)
    double chi_grad = i_generator->fChiQuadro(par);
    double chi_prec = i_generator->fChiQuadro(par_prec);
    if (chi_grad < chi_prec)
    {
        //par sono già modificati e i migliori
        chi_quadro_min = chi_grad;
    }
    else
    {
        par = par_prec;                 //Riporto i parametri e il chi_quadro ai valori originali (prima della chiamata della funzione)
        chi_quadro_min = chi_prec;      //
    }

    return;
}


// gradiente della funzione chi quadro
std::vector<double> gradient_descent_algorithm::grad_f_chi_quadro(const std::vector<double> par) {
    std::vector<double> grad;
    int par_size = par.size();

    // Calcolo della derivata parziale per ogni parametro
    for (int i = 0; i < par_size; ++i) {
        
        double h = 1e-8 * (fabs(par[i]) < 1) ? 1 : fabs(par[i]);      // 1e-8 è circa sqrt(epsilon) con epsilon= 2.22e-16 che è la precisione in double

        // Crea una copia del vettore dei parametri
        std::vector<double> par_plus = par;
        std::vector<double> par_minus = par;

        // Modifica il parametro i-esimo in positivo e negativo
        par_plus[i] += h;
        par_minus[i] -= h;

        // Approssima la derivata con la differenza centrale
        grad.emplace_back((i_generator->fChiQuadro(par_plus) - i_generator->fChiQuadro(par_minus)) / (2 * h));
        //cout <<par_plus[i] << "\t" << par_minus[i] << "\t" << i_generator->fChiQuadro(par_plus) << "\t" << i_generator->fChiQuadro(par_minus) << "\t" << grad[i] << endl;
    }

    return grad;
}


gradient_descent_algorithm::~gradient_descent_algorithm() {
    delete i_generator;
}