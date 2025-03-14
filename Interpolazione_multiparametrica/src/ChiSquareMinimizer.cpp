#include "ChiSquareMinimizer.h"

#include "interpolator.h"

#include "covering.h"
#include "results.h"

#include "matrix.h"
#include "chi_square.h"
#include "gradient_descent_algorithm.h"
#include "util.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <thread>
#include <chrono>

ChiSquareMinimizer::ChiSquareMinimizer(std::map<std::string, bool>& options):
    faster(options["faster"]), approx(options["approx"]), complex(options["complex"]), retta(options["retta"])
{

}

void ChiSquareMinimizer::begin(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    double sensibility = 0.1;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);

    if (faster || complex) {        //Stampa subito a schermo dei risultati (se in modalità 'complex' oppure 'faster')
        Results(par_best, approx, std::cout);       // passando 'std::cout' stampo a schermo, altrimenti potrei passargli un ofstream perché scriva su un file ad esempio
    }
}


void ChiSquareMinimizer::compute(std::vector<double>& par_best) {

    std::vector<double> par_best_prec(par_best.size(), 0.0);


    for (int cicle_programms = 1; cicle_programms < 10; cicle_programms++)
    {

        Interpolator* i_generator = Interpolator::getNewInstance();        // deve riaggiornarsi qui dentro ogni volta

        computeCovering(par_best, cicle_programms);

        double sensibility1 = 0.01;
        gradient_descent_algorithm(par_best, chi_quadro_min, sensibility1);

        /*
        std::cout << "Ora chiamo Results(" << std::endl;
        for (int i = 0; i < par_best.size(); i++)
        {
            std::cout << "par" << i << " " << par_best[i] << "\t";
        }
        std::cout << std::endl;
        */

        if (i_generator->fChiQuadro(par_best_prec) - i_generator->fChiQuadro(par_best) < 0.001) {
            delete i_generator;
            return;
        }
        else if (complex) {     // se non sono alla fine e voglio risultati intermedi ('complex=true'), allora stampo i risultati a schermo
            Results(par_best, approx, std::cout);
        }
        delete i_generator;

        par_best_prec = par_best;
    }
    
}


void ChiSquareMinimizer::end(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente' prima di dare i risultati definitivi
    double sensibility1 = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility1);

    // stampo i risultati a schermo
    Results(par_best, approx, std::cout);

}




void ChiSquareMinimizer::computeCovering(std::vector<double>& par_best, int cicle_programms) {
    
    c_generator = new covering(par_best, chi_quadro_min, cicle_programms, complex, faster);


    int cicle = 1;      //primo ciclo di livelli
    for (int k = 0; k < 10; k++) {

        //Stamapa a schermo a che ciclo di programma si è arrivati
        c_generator->status(cicle, k);

        // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
        c_generator->ricoprimento(0, false);   //I migliori parametri sono in 'par_best'

        // Passaggio al livello successivo
        c_generator->next();
    }


    int num_level_tasks = 4;        // faccio il livello 0 (singolo cubo n-dim attorno al punto), livello 1 (primo ricoprimento attorno al primo cubo n-dim), livello 2, ecc. fino al livello n = 'num_level_tasks'
    unsigned int min_tasksCompleted = std::pow(2 * num_level_tasks + 1, par_best.size());
    unsigned int time_max = 15000;

    if (cicle_programms > 1) {
        num_level_tasks = 4;
        time_max = 5000;
    }

    // Devo aspettare che il pool di multi-thread abbia completato tutte le tasks
    c_generator->GetResults(par_best, 1000, min_tasksCompleted, time_max); // in ordine: parametri su cui salvare i migliori, tempo minimo analisi (es. 1s), minimo numero tasks (es. (2*4+1)^3=729), tempo massimo analisi (es. 15s)
    
    delete c_generator;
    
}