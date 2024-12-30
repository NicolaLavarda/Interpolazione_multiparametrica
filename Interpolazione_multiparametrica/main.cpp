#include "file.h"
#include "input.h"
#include "automatic_research.h"
#include "covering.h"
#include "linear_mode.h"
#include "matrix.h"
#include "interpolating_function.h"
#include "results.h"
#include "plot.h"
#include "chi_square.h"
#include "gradient_descent_algorithm.h"
#include "util.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <exception>
#include <functional>
#include <map>

#include "TCanvas.h"
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>

using namespace std;

//Variabili globali
vector<double> x, sigma_x, y, sigma_y;                    //Dati iniziali
double chi_quadro_min = 1e30;                             //Chi quadro minimo assoluto (in ogni istante di tutto il programma)
vector<double> chi_quadro;                                //Vettore chi quadri minimi (ripulito ad ogni riesecuzione del programma)
vector<double> par_best;                                  //Parametri effettivamente stampati poi a schermo
vector<int> livelli(100, 0);                              //Livelli di ogni ricoprimento (valori diversi da 0 se in quel livello viene migliorato il chi quadro)
vector<vector<double>> par_matrix;                        //par_best di ogni ciclo di riesecuzione del programma
int cicle_programms = 0;                                  //Numero di cicli di riesecuzione del programma


int main(int argc, char* argv[]) {

    //------------LETTURA INPUT-------------------------

    vector<double> par;
    int num_a = 0;

    std::map<std::string, bool> options = {
        {"faster", false},
        {"approx", false},
        {"complex", false},
        {"retta", false},
        {"save", false},
        {"plot", false}
    };

    input(argc, argv, par, num_a, options);

    bool faster = options["faster"] ? true : false;
    bool approx = options["approx"] ? true : false;
    bool complex = options["complex"] ? true : false;
    bool retta = options["retta"] ? true : false;
    bool save = options["save"] ? true : false;
    bool plot = options["plot"] ? true : false;

    
    //------------RICERCA AUTOMATICA--------------------

    if (num_a != 0)
        parametri_auto(par, chi_quadro_min, complex);
    par_best = par;


    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    double sensibility = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);

    if (faster || complex) {        //Stampa subito a schermo dei risultati (se in modalità 'complex' oppure 'faster')
        Results(par_best, approx, std::cout);       // passando 'std::cout' stampo a schermo, altrimenti potrei passargli un ofstream perché scriva su un file ad esempio
    }


    //------------ESECUZIONE PROGRAMMA------------------

    //eseguo il blocco seguente fino ad una effettiva variazione del chi quadro minore del per mille (0.001) -> condizione presente alla fine del ciclo while in 'end()'
    
    bool errore_lin = false;

    for (int cicle_programms = 1; cicle_programms < 20; cicle_programms++) {  //la condizione di termine è in 'end()' [...] '<20' è giusto per essere sicuro che per qualche bug non vada all'infinito
        
        

        covering c_generator(par_best, chi_quadro_min, cicle_programms, complex, faster);
        linear_mode l_generator(par_best, faster, complex);

        //Calcolo parametri
        int cicle = 1;      //primo ciclo di livelli
        for (int k = 0; k < 15; k++) {

            //Stamapa a schermo a che ciclo di programma si è arrivati
            c_generator.status(cicle, k);
            
            // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
            c_generator.ricoprimento(par, par_best, 0, false);   //I migliori parametri sono in 'par_best'


            if (complex)    // l'analisi è presumibilmente lunga e difficile (forse instabile) quindi voglio vedere i parametri migliorati ad ogni fine di livello [-> volendo, aggiungere all'if '&& cicle_programms == 1']
            {
                cout << "\t -> ";
                for (int k = 0; k < par.size(); k++)
                    cout << par_best[k] << "\t";
            }


            //Modalità con ricerca lungo una retta di miglioramento dei parametri (utile se sono molto distante)
            bool ricerca_retta = false;     //Cambia in 'true' se il metodo qui sotto della ricerca lungo la retta funziona
            if (par.size() == 1) errore_lin = true;
            if ((!errore_lin || retta) && cicle == 1)       // se da input metto retta=true vuol dire che voglio che venga sempre usato quando possibile il metodo della retta
                l_generator.research(par_best, errore_lin, ricerca_retta);      //fin tanto che non si hanno almeno un tot di punti prefissati li si raccolgono, poi 'linear_mode' continua effettivamente cercando di migliorare i parametri
            

            // Passaggio al livello o al ciclo successivo
            c_generator.next();

            // Motivi di uscita dal ciclo 'for'
            if (c_generator.exit(ricerca_retta)) break;

        }

        // Miglioro i parametri usando il metodo 'discesa_gradiente'
        double sensibility = 0.01;
        gradient_descent_algorithm(par_best, chi_quadro_min, sensibility);

        // Se non esco allora stampo i risultati (se 'complex' è attivato)
        if (complex)
        {
            //Stampo a schermo i risultati
            Results(par_best, approx, std::cout);
        }

        // riaggiorno parametri, libero 'chi_quadro' e 'livelli', e passo al ciclo di programma successivo [...] poi in caso esco e termino tutto
        if (c_generator.end()) break;

    }

    // Miglioro i parametri usando il metodo 'discesa_gradiente'
    double sensibility1 = 0.01;
    gradient_descent_algorithm(par_best, chi_quadro_min, sensibility1);

    // stampo i risultati a schermo una sola volta alla fine se 'complex=false' altrimenti sono già stati stampati a schermo
    if (!complex)
    {
        //Stampo a schermo i risultati
        Results(par_best, approx, std::cout);
    }

    //------------SALVATAGGIO DEI RISULTATI-----------------

    if (save)
    {
        //salvo i risultati nel file di testo
        writeFile(string(argv[1]), x, sigma_x, y, sigma_y, par_best, approx, string(argv[par_best.size() + 2]));
    }


    //------------PRODUZIONE DEI GRAFICI--------------------

    if (plot)
    {
        PlotGenerator generator(par_best, x, sigma_x, y, sigma_y, string(argv[1]));
        generator.compute_plot_function();
        generator.compute_plot_chi_distribution();
        generator.compute_plot_Residuals();
        generator.save(".png");
    }


    return 0;
}