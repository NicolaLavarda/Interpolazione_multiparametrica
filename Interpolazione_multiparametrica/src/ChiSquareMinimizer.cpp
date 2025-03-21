#include "ChiSquareMinimizer.h"

//#include "interpolator.h"

#include "matrix.h"
#include "chi_square.h"
#include "gradient_descent_algorithm.h"
#include "BranchingCover.h"
#include "covering.h"
#include "results.h"
#include "util.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <thread>
#include <chrono>

ChiSquareMinimizer::ChiSquareMinimizer(std::map<std::string, bool>& options):
    faster(options["faster"]), approx(options["approx"]), complex(options["complex"])
{

}

void ChiSquareMinimizer::begin(std::vector<double>& par_best) {
        
    // Miglioro i parametri usando il metodo 'BranchingCover'
    BranchingCover b_generator2(par_best);
    b_generator2.SetModel(7, 1, 0.05);
    b_generator2.compute(par_best);

    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    gradient_descent_algorithm(par_best, chi_quadro_min, 0.01);

    // Normalizzo i parametri perché abbiano lo stesso ordine di grandezza
    Interpolator::normalizeTo10(par_best);

    if (faster || complex) {        //Stampa subito a schermo dei risultati (se in modalità 'complex' oppure 'faster')
        Results(par_best, approx, std::cout);       // passando 'std::cout' stampo a schermo, altrimenti potrei passargli un ofstream perché scriva su un file ad esempio
    }
}


void ChiSquareMinimizer::compute(std::vector<double>& par_best) {

    std::vector<double> par_best_prec(par_best.size(), 0.0);
    double chi_quadro_min_prec = chi_quadro_min * 100;      // dev'essere maggiore all'inizio ("simula un chi quadro precedente iniziale")

    for (int cicle_programms = 1; cicle_programms < 10; cicle_programms++)
    {
        // Normalizzo i parametri perché abbiano lo stesso ordine di grandezza
        Interpolator::normalizeTo10(par_best);

        // miglioro i parametri con un "ricoprimento" dello spazio dei parametri (vedi 'covering')
        computeCovering(par_best, cicle_programms);

        // miglioro i parametri con un algoritmo di discesa del gradiente
        gradient_descent_algorithm(par_best, chi_quadro_min, 0.005);         // sensibility = 0.005

        // se sono migliorati entro solo 5 sigma allora diminuisco i passi di 'covering' a 1 sigma
        decreseSteps(par_best, par_best_prec);

        if (chi_quadro_min_prec - chi_quadro_min < 0.001) {
            return;
        }
        else if (complex) {     // se non sono alla fine e voglio risultati intermedi ('complex=true'), allora stampo i risultati a schermo
            Results(par_best, approx, std::cout);
        }

        par_best_prec = par_best;
        chi_quadro_min_prec = chi_quadro_min;
    }
    
}


void ChiSquareMinimizer::end(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente' prima di dare i risultati definitivi
    gradient_descent_algorithm(par_best, chi_quadro_min, 0.001);

    // stampo i risultati a schermo
    Results(par_best, approx, std::cout);

}


void ChiSquareMinimizer::decreseSteps(std::vector<double>& par_best, std::vector<double>& par_best_prec) {

    // Controllo se i parametri attuali differiscono da quelli del ciclo precedente entro 3 sigma dei presenti parametri
    std::vector<double> errors = Results::GetErrors(par_best);
    int check = 0;
    int n = par_best.size();
    for (int i = 0; i < n; i++)
    {
        if (std::fabs(par_best[i] - par_best_prec[i]) < 5 * errors[i])
            check++;
    }
    if (check == n)     // i parametri non sono migliorati più di 3 sigma dal ciclo precedente
    {
        // Riduco il passo del prossimo covering su ogni parametro pari a 1 sigma dell'errore sul rispettivo parametro
        steps = errors;     // poi utilizzo 'c_generator->SetSteps(steps)' in 'computeCovering'
        update_steps = true;
    }
    else
        update_steps = false;
}


void ChiSquareMinimizer::computeCovering(std::vector<double>& par_best, int cicle_programm) {
    
    c_generator = new covering(par_best, chi_quadro_min, cicle_programm, complex);

    // se sono nelle condizioni definite in 'compute' definisco dei passi personalizzati per il ricoprimento
    if (update_steps)
        c_generator->SetSteps(steps);

    //Stamapa a schermo a che ciclo di programma si è arrivati
    c_generator->status();

    int num_level_tasks = 4;        // faccio il livello 0 (singolo cubo n-dim attorno al punto), livello 1 (primo ricoprimento attorno al primo cubo n-dim), livello 2, ecc. fino al livello n = 'num_level_tasks'
    for (int k = 0; k < num_level_tasks + 1; k++) {

        // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
        c_generator->ricoprimento(0, false);   //I migliori parametri sono in 'par_best'

        // Passaggio al livello successivo
        c_generator->next();
    }

    unsigned int min_tasksCompleted = std::pow(2 * num_level_tasks + 1, par_best.size());
    unsigned int time_max = 5000;

    if (cicle_programm > 1)
        time_max = 2000;

    // Devo aspettare che il pool di multi-thread abbia completato tutte le tasks
    c_generator->GetResults(par_best, 1000, min_tasksCompleted, time_max); // in ordine: parametri su cui salvare i migliori, tempo minimo analisi (es. 1s), minimo numero tasks (es. (2*4+1)^3=729), tempo massimo analisi (es. 15s)
    
    delete c_generator;
    
}