#ifndef BISECTION_ALGORITHM_H
#define BISECTION_ALGORITHM_H

#include "interpolator.h"

#include <vector>
#include <memory>

// Chiamare nel programma 'bisection_algorithm(par, step, chi_quadro_min);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


class bisection_algorithm {
public:

    bisection_algorithm(std::vector<double>& par, const std::vector<double> step, double& chi_quadro_min);

    ~bisection_algorithm();

private:

    void bisection(std::vector<double> par, std::vector<double>& par_def, const std::vector<double> passo, int n);

    // Accesso ad una nuova istanza già settata (per non rallentare i thread aspettandosi a vicenda i mutex lock per accedere uno alla volta a 'fChiQuadro')
    Interpolator* i_generator = Interpolator::getNewInstance();
};



#endif