#include <file.h>
#include <input.h>
#include <interpolating_function.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <regex>
#include <functional>
#include <map>

using namespace std;

extern vector<double> x, sigma_x, y, sigma_y;       //Dati iniziali
extern bool flag;   //considerare o meno anche errori in x

bool isNumber(const string& str) {
    // Regex per un numero valido (intero o floating point, anche esponenziale)
    // -> Numeri validi: 
    //      Numero intero (123, -456)
    //      Numero decimale (123.45, -0.67, .89)
    //      Notazione scientifica (1.23e4, -4.56E-7)

    regex numberPattern(R"(^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$)");
    return regex_match(str, numberPattern);
}

void input(int argc, char* argv[], vector<double>& par, int& num_a, map<string, bool>& options) {


    if (argc == 1) {
        cerr << "No file name has been entered" << endl;
        exit(EXIT_FAILURE);
    }
    int num_columns = readFile(string(argv[1]), x, sigma_x, y, sigma_y);
    if (num_columns == 4) flag = false;     //considero errori anche sulle x
    if (num_columns < 0)
        exit(EXIT_FAILURE);

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
        setup_expression(espressione_interpolante);
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

