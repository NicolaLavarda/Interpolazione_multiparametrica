#include "input.h"

#include "util.h"
#include "file.h"
#include "interpolator.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <functional>
#include <map>

using namespace std;


void input(int argc, char* argv[], std::vector<double>& par,
    int& num_a, std::map<std::string, bool>& options,
    std::vector<double>& x, std::vector<double>& sigma_x,
    std::vector<double>& y, std::vector<double>& sigma_y) {
    
    // Prendo l'istanza al Singleton 'Interpolator'
    Interpolator& i_generator = Interpolator::getInstance();

    if (argc == 1) {
        cerr << "No file name has been entered" << endl;
        exit(EXIT_FAILURE);
    }
    
    //std::vector<double> x, sigma_x, y, sigma_y;
    int num_columns = readFile(string(argv[1]), x, sigma_x, y, sigma_y);
    if (num_columns < 0)
        exit(EXIT_FAILURE);
    i_generator.setData(x, sigma_x, y, sigma_y);

    if (argc == 2) {
        cerr << "No parameters value have been entered" << endl;
        exit(EXIT_FAILURE);
    }


    //Inizializzazione vettore parametri (ricerca logaritmica automatica o valori attorno cui cercare inseriti in compilazione dall'utente)
    // ricordo che 'vector<double> par;' e 'int num_a = 0;' sono definiti nella main e passati qui come argomenti
    
    for (int i = 2; i < argc; i++) {
        try {
            if (isNumber(string(argv[i])))
                par.push_back(stod(argv[i]));
            else
                throw invalid_argument("Non convertibile in numero");
        }
        catch (const invalid_argument& e) {
            if (string(argv[i]) == "a") {
                par.push_back(0);           //per far capire alle funzioni di ricerca automatica che i parametri pari esattamente a 0 sono da ricercare automatici
                num_a++;
            }
            else
                break;          //if (string(argv[i]) != "a") break;
        }
    }
    int par_size = par.size();


    if (argc == par_size + 2) {
        cerr << "No interpolating function has been entered" << endl;
        exit(EXIT_FAILURE);
    }


    //Inizializzazione e interpretazione della funzione interpolante (dev'essere inserita da terminale tra virgolette "...")
    string espressione_interpolante = string("(1)*") + argv[par_size + 2];        //"(1)*" è perché se scrivo come primo carattere un numero mi da errore e non capisco perché (in questo modo ho risolto)
    try {
        // Inizializza l'espressione solo una volta
        i_generator.setExpression(espressione_interpolante);
    }
    catch (const runtime_error& e) {
        std::cerr << "Error in compilation of the expression: " << espressione_interpolante << std::endl;
        cerr << "Error: " << e.what() << endl;
    }


    //Parametri ulteriori inseriti come ultima cosa in comando di compilazione dall'utente
    if (argc > par_size + 3)
    {
        for (int i = par_size + 3; i < argc; ++i) {
            std::string arg = argv[i];
            // Controlla se l'argomento è una chiave valida nella mappa
            if (options.find(arg) != options.end()) {
                options[arg] = true; // Imposta il valore associato a true
            }
            else {
                std::cerr << "Not recognized flag: " << arg << endl;
            }
        }
    }


}

