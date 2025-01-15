#ifndef BISECTION_ALGORITHM_H
#define BISECTION_ALGORITHM_H

#include "interpolator.h"

#include <vector>

// Chiamare nel programma 'bisection_algorithm(par, step, chi_quadro_min);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


class bisection_algorithm {
public:

    bisection_algorithm(std::vector<double>& par, const std::vector<double> step, double& chi_quadro_min);

private:

    void bisection(std::vector<double> par, std::vector<double>& par_def, const std::vector<double> passo, int n);

    // Riferimento all'istanza Singleton
    Interpolator& i_generator = Interpolator::getInstance();
};



#endif