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

void ChiSquareMinimizer::compute(std::vector<double>& par_best) {

    //eseguo il blocco seguente fino ad una effettiva variazione del chi quadro minore del per mille (0.001) -> condizione presente alla fine del ciclo while in 'end()'

    for (int cicle_programms = 1; cicle_programms < 20; cicle_programms++) {  //la condizione di termine è in 'end()' [...] '<20' è giusto per essere sicuro che per qualche bug non vada all'infinito

        covering c_generator(par_best, chi_quadro_min, cicle_programms, complex, faster);
        linear_mode l_generator(par_best, faster, complex);

        //Calcolo parametri
        int cicle = 1;      //primo ciclo di livelli
        for (int k = 0; k < 15; k++) {

            //Stamapa a schermo a che ciclo di programma si è arrivati
            c_generator.status(cicle, k);

            // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
            c_generator.ricoprimento(par_best, 0, false);   //I migliori parametri sono in 'par_best'


            if (complex)    // l'analisi è presumibilmente lunga e difficile (forse instabile) quindi voglio vedere i parametri migliorati ad ogni fine di livello [-> volendo, aggiungere all'if '&& cicle_programms == 1']
            {
                cout << "\t -> ";
                for (int k = 0; k < par_best.size(); k++)
                    cout << par_best[k] << "\t";
            }

            //Modalità con ricerca lungo una retta di miglioramento dei parametri (utile se sono molto distante)
            bool ricerca_retta = false;     //Cambia in 'true' se il metodo qui sotto della ricerca lungo la retta funziona
            if (par_best.size() == 1) errore_lin = true;    //non ha senso se c'è un solo parametro
            if ((!errore_lin || retta) && (!faster && cicle == 1))       // se da input metto 'retta=true' vuol dire che voglio che venga sempre usato quando possibile il metodo della retta, con 'faster=true' l'analisi è verosimilmente facile e non perdo tempo con la retta
                l_generator.research(par_best, errore_lin, ricerca_retta);      //fin tanto che non si hanno almeno un tot di punti prefissati li si raccolgono, poi 'linear_mode' continua effettivamente cercando di migliorare i parametri


            // Passaggio al livello o al ciclo successivo
            c_generator.next();

            // Motivi di uscita dal presente ciclo 'for'
            if (c_generator.exit(ricerca_retta)) break;

        }

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
        if (c_generator.end()) break;

    }

}

void ChiSquareMinimizer::end(std::vector<double>& par_best) {

    // Miglioro i parametri usando il metodo 'discesa_gradiente' prima di dare i risultati definitivi
    double sensibility1 = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility1);

    // stampo i risultati a schermo una sola volta alla fine se 'complex=false' altrimenti sono già stati stampati a schermo
    if (!complex)
        Results(par_best, approx, x_size, std::cout);

}

