#include "BranchingCover.h"

#include "interpolator.h"
#include "gradient_descent_algorithm.h"

#include <iostream>
#include <vector>
#include <algorithm>  // per std::max_element
#include <cmath>

BranchingCover::BranchingCover(std::vector<double>& par_best):
    par_size(par_best.size()), level(0), max_level(5)
{
    // Mi assicuro di partire da un minimo locale all'inizio
    double sensibility = 0.1;
    double chi_quadro_min = 0.0;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);
}

void BranchingCover::doBetter(std::vector<double>& par_best) {
    doBetter(par_best, std::pow(10, max_level));
}

void BranchingCover::doBetter(std::vector<double>& par_best, int address) {
    level++;

    // creo 'div' insieme di gradienti che puntano in diverse direzioni
    std::vector<std::vector<double>> div = division_direction(par_best);

    int n = div.size();
    for (int i = 0; i < n; i++)
    {
        std::vector<double> par_i = par_best;       // 'par_best' punto da cui partire
        up_down(par_i ,div[i]);    // sale finché non trova una "buca" e lì allora scende finio al minimo locale
        par_div.push_back(par_i);   // arrivato in un minimo locale lo salvo
        int ad = register_par_address(i, address);      // associo alla posizione del vettore nella ramificazione l'indice corrispondente in 'par_div' (in questo momento è oviamente l'ultimo essendo appena stato aggiunto)
        if (level != max_level)
            doBetter(par_i, ad);   // si divide in più direzioni e riparte
    }
    level--;
}

int BranchingCover::register_par_address(int i, int address) {
    // 'address' che riceve è quello dei 'par' da cui è stato generato
    address += std::pow(10, max_level - level) * (i + 1);       // 12140 significa che è nel 4° ramo del 1° ramo del 2° ramo del 1° ramo (e quindi ad una profondità 4 su 5 che è quella massima)
        
    address_par[address] = par_div.size() - 1;
    std::cout << "address: " << address << std::endl;
    return address;
}


bool BranchingCover::up_down(std::vector<double>& par_i, std::vector<double> grad) {

    std::vector<double> par_prov = par_i;
    
    double chi_min = i_generator->fChiQuadro(par_i);
    double chi_min0 = chi_min;
    double tasso_apprendimento = 0.01;



    std::vector<double> check_max;
    for (int i = 0; i < par_size; i++)
    {
        check_max.push_back(fabs(grad[i] / par_i[i]));
        //cout << " -> " << fabs(grad[i] / par[i]);
    }
    //cout << endl;
    auto max_iter = std::max_element(check_max.begin(), check_max.end());
    int max_index = std::distance(check_max.begin(), max_iter);

    /*
    // trovo l'indice del valore massimo tra i 'grad[i]/par[i]'
    auto max_iter = std::max_element(grad.begin(), grad.end(),
        [&](double a, double b) { return fabs(a / par_i[&a - &grad[0]]) < fabs(b / par_i[&b - &grad[0]]); });

    int max_index = max_iter - grad.begin();  // Ottiene direttamente l'indice (come distanza tra gli iteratori)
    */
    //std::cout << max_index << std::endl;

    for (int iter = 0; iter < 100; ++iter) {

        std::vector<double> par_p = par_i;
        for (int i = 0; i < par_size; ++i) {
            par_p[i] += grad[i] * par_p[i] * tasso_apprendimento / fabs(grad[max_index]);     // il "+=" è perché voglio risalire (controgradiente per uscire dal minimo locale)
        }

        double chi_min_p = i_generator->fChiQuadro(par_p);

        std::cout << "Iterazione " << iter + 1 << ": i_generator->fChiQuadro(par) = " << i_generator->fChiQuadro(par_p); for (int i = 0; i < par_size; i++) std::cout << "\t" << par_p[i]; std::cout << std::endl;

        if (chi_min_p > chi_min0)
        {
            if (chi_min_p / chi_min0 > 1.02)    // se il chi_quadro aumenta più del 2% lo rallento dimezzando 'tasso_apprendimento' (rischia di salire a salti troppo grandi e non mi becca più "buche" di minimi su cui ricadere)
            {
                tasso_apprendimento /= 2;
            }
            chi_min0 = chi_min_p;
            par_prov = par_p;
        }
        else   // sono caduto in una buca in cui cercare un minimo
        {
            chi_min0 = chi_min_p;
            par_prov = par_p;
            break;  // termino la salita per iniziare la discesa
        }
        //std::cout << "Iterazione " << iter + 1 << ": i_generator->fChiQuadro(par) = " << i_generator->fChiQuadro(par_prov); for (int i = 0; i < par_size; i++) std::cout << "\t" << par_prov[i]; std::cout << std::endl;
    }

    gradient_descent_algorithm(par_prov, chi_min0, tasso_apprendimento);     // vado ora in discesa

    std::cout << "Final: i_generator->fChiQuadro(par) = " << i_generator->fChiQuadro(par_prov); for (int i = 0; i < par_size; i++) std::cout << "\t" << par_prov[i]; std::cout << std::endl;

    return (chi_min0 < chi_min) ? true : false;     // 'true' se ha migliorato il chi_quadro

}

// ritorna n gradienti differenti per andare in direzioni diverse
std::vector<std::vector<double>> BranchingCover::division_direction(const std::vector<double> par) {
    std::vector<std::vector<double>> div;

    int n = par.size();
    for (int i = 0; i < n; i++) {
        std::vector<double> div_n = grad_f_chi_quadro(par);     // Il gradiente rimane costante per continuare a salire sempre in quella direzione (finché non comincia a scendere in una buca di minimo)
        div_n[i] = 0;       // cambiano solo gli altri n-1 parametri (divido in n direzioni)
        div.emplace_back(div_n);
    }

    return div;
}






// gradiente della funzione chi quadro
std::vector<double> BranchingCover::grad_f_chi_quadro(const std::vector<double> par) {
    std::vector<double> grad;
    int par_size = par.size();

    // Calcolo della derivata parziale per ogni parametro
    for (int i = 0; i < par_size; ++i) {

        double h = 1e-8 * fabs(par[i]);      // 1e-8 è circa sqrt(epsilon) con epsilon= 2.22e-16 che è la precisione in double

        // Crea una copia del vettore dei parametri
        std::vector<double> par_plus = par;
        std::vector<double> par_minus = par;

        // Modifica il parametro i-esimo in positivo e negativo
        par_plus[i] += h;
        par_minus[i] -= h;

        // Approssima la derivata con la differenza centrale
        grad.emplace_back((i_generator->fChiQuadro(par_plus) - i_generator->fChiQuadro(par_minus)) / (2 * h));
    }

    return grad;
}


BranchingCover::~BranchingCover() {
    delete i_generator;
}