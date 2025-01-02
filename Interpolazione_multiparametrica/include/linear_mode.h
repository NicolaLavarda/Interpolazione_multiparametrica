#ifndef LINEAR_MODE_H
#define LINEAR_MODE_H

#include "interpolator.h"

#include <vector>

// Chiamare nel programma 'linear_mode(par_lin, m_lin, q_lin, errore_lin, ricerca_retta, faster, complex);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


class linear_mode {
public:

    //costruttore di base
    linear_mode(std::vector<double> par, bool faster, bool complex);

    //funzione effettiva da usare nella main
    void research(std::vector<double> par_in, bool& errore_lin, bool& ricerca_retta);

private:

    void linearFit(std::vector<double>& x, std::vector<double>& y, double& m, double& q);

    void bisezione_lin(std::vector<double>& par, std::vector<double> m, std::vector<double> q, int n_worst);

    // Funzione per trovare il parametro "peggiore" -> devo capire qual è il parametro peggiore in funzione del quale migliorare gli altri con 'linear_mode'
    int worst_parameter(std::vector<double>& m_lin, std::vector<double>& q_lin);


    double chi_quadro_min;              //Chi quadro minimo assoluto
    std::vector<double> par_best;       //Attuali migliori parametri

    //Blocco inizializzazioni per ricerca lungo retta -> metodo della retta (molto più efficente dei "ricoprimenti" quando sono distante dalla soluzione)
    std::vector<std::vector<double>> par_lin;              //parametri da interpolare linearmente
    std::vector<double> m_lin;     //coefficienti angolari per parametri da interpolare linearmente
    std::vector<double> q_lin;     //intercette per parametri da interpolare linearmente

    bool faster;
    bool complex;

    // Riferimento all'istanza Singleton
    Interpolator& i_generator = Interpolator::getInstance();

};

#endif