#include "ChiSquareMinimizer.h"

#include "covering.h"
#include "linear_mode.h"
#include "results.h"

#include "matrix.h"
#include "chi_square.h"
#include "gradient_descent_algorithm.h"
#include "util.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include <thread>
#include <chrono>

ChiSquareMinimizer::ChiSquareMinimizer(std::map<std::string, bool>& options, int num_val):
    faster(options["faster"]), approx(options["approx"]), complex(options["complex"]), retta(options["retta"]), x_size(num_val)
{

}

void ChiSquareMinimizer::begin(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    double sensibility = 0.1;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);

    if (faster || complex) {        //Stampa subito a schermo dei risultati (se in modalità 'complex' oppure 'faster')
        Results(par_best, approx, x_size, std::cout);       // passando 'std::cout' stampo a schermo, altrimenti potrei passargli un ofstream perché scriva su un file ad esempio
    }

}



void ChiSquareMinimizer::thread1(std::vector<double>& par_best) {

    std::cout << "Attesa di 1000 millisecondi..." << std::endl;

    // Aspetta 1000 millisecondi
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    std::cout << "Attesa completata!" << std::endl;

}


void ChiSquareMinimizer::thread2(std::vector<double>& par_best) {

    int cicle_programms = 1;
    c_generator = new covering(par_best, chi_quadro_min, cicle_programms, complex, faster);


    //Calcolo parametri
    int cicle = 1;      //primo ciclo di livelli
    for (int k = 0; k < 10; k++) {

        //Stamapa a schermo a che ciclo di programma si è arrivati
        c_generator->status(cicle, k);

        // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
        c_generator->ricoprimento(par_best, 0, false);   //I migliori parametri sono in 'par_best'


        /*
        if (complex)    // l'analisi è presumibilmente lunga e difficile (forse instabile) quindi voglio vedere i parametri migliorati ad ogni fine di livello [-> volendo, aggiungere all'if '&& cicle_programms == 1']
        {
            cout << "\t -> ";
            for (int k = 0; k < par_best.size(); k++)
                cout << par_best[k] << "\t";
        }
        */

        // Passaggio al livello o al ciclo successivo
        c_generator->next();


        //std::cout << "finishhhhhhhhhhh = " << c_generator->GetFinish() << std::endl;

        
        /*
        if (c_generator->GetFinish())
        {
            std::cout << "escooooo" << std::endl;
            return;
        }
        */
        
    }
}



void ChiSquareMinimizer::compute(std::vector<double>& par_best) {

    /*
    // Creazione dei due thread
    //std::thread t1([this, &par_best]() { this->thread1(par_best); });
    std::thread t2([this, &par_best]() { this->thread2(par_best); });

    std::cout << "Attesa di 4 secondi..." << std::endl;

    // Aspetta 4000 millisecondi
    std::this_thread::sleep_for(std::chrono::milliseconds(4000));

    std::cout << "Attesa completata!" << std::endl;

    c_generator->SetFinish(true);

    // Aspetta la terminazione del thread
    t2.join();
    */

    thread2(par_best);


    //std::this_thread::sleep_for(std::chrono::milliseconds(4000));

    // Devo aspettare che il thread finisca. In particolare vuol dire che anche il pool di multi-thread ha completato tutte le tasks
    c_generator->GetResults(par_best, 1000, 50);


    /*
    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    double sensibility = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);

    // Se non esco allora stampo i risultati (se 'complex' è attivato)
    if (complex)
    {
        //Stampo a schermo i risultati
        Results(par_best, approx, x_size, std::cout);
    }

    // riaggiorno parametri, libero 'chi_quadro' e 'livelli', e passo al ciclo di programma successivo [...] poi in caso esco e termino tutto
    if (c_generator->end()) break;
    */
}



void ChiSquareMinimizer::end(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente' prima di dare i risultati definitivi
    double sensibility1 = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility1);

    // stampo i risultati a schermo
    Results(par_best, approx, x_size, std::cout);

}

