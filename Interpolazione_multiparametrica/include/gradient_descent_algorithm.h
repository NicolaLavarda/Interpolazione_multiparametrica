#ifndef GRADIENT_DESCENT_ALGORITHM_H
#define GRADIENT_DESCENT_ALGORITHM_H

#include "interpolator.h"

#include <vector>
#include <memory>

// Chiamare nel programma 'gradient_descent_algorithm(par_best, approx);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


// Classe Base
class gradient_descent_algorithm {
public:

    gradient_descent_algorithm(std::vector<double>& par, double& chi_quadro_min, double sensibility);

    ~gradient_descent_algorithm();

private:

    std::vector<double> grad_f_chi_quadro(const std::vector<double> par);

    // Accesso ad una nuova istanza già settata
    Interpolator* i_generator = Interpolator::getNewInstance();

};



#endif