#include "input.h"

#include "util.h"
#include "file.h"
#include "interpolator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <functional>
#include <map>

using namespace std;

input::input(int argc, char* argv[])
    : argc(argc), argv(argv) {}


void input::compute(std::vector<double>& par, int& num_a, std::map<std::string, bool>& options,
                     std::vector<double>& x, std::vector<double>& sigma_x,
                     std::vector<double>& y, std::vector<double>& sigma_y) {

    if (argc == 1) {
        cerr << "No file name has been entered" << endl;
        exit(EXIT_FAILURE);
    }

    filePath = string(argv[1]);
    int num_columns = readFile(filePath, x, sigma_x, y, sigma_y);
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
    for (auto& [key, value] : options) {
        value = false;
    }
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


bool input::improved(const std::string filePath, std::vector<double> par_best) {
    return (GetValueFromFile(filePath, "chi-squared", "=") > i_generator.fChiQuadro(par_best)) ? true : false;
}


std::vector<double> input::GetParametersFromFile() {
    std::vector<double> value_parameters;
    for (size_t i = 0; i < name_parameters.size(); i++)
    {
        try
        {
            value_parameters.emplace_back(GetParameter(name_parameters[i]));
        }
        catch (const std::runtime_error& e)
        {
            break;
        }
    }
    return value_parameters;
}


std::vector<double> input::GetErrorsFromFile() {
    std::vector<double> value_errors;
    for (size_t i = 0; i < name_parameters.size(); i++)
    {
        try
        {
            value_errors.emplace_back(GetErrorOfParameter(name_parameters[i]));
        }
        catch (const std::runtime_error& e)
        {
            break;
        }
    }
    return value_errors;
}


double input::GetParameter(const std::string parameter) {
    return GetValueFromFile(filePath, parameter + " =", "(");
}


double input::GetErrorOfParameter(const std::string parameter) {
    return GetValueFromFile(filePath, parameter + " =", "+-");
}



double input::GetValueFromFile(const std::string filePath, const std::string target, const std::string afterThat) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file: " << filePath << std::endl;
        throw runtime_error("Unable to open the file");
        return -1; // codice di errore per file non trovato
    }

    std::string line;
    while (std::getline(file, line)) {
        // Cerca la riga che contiene 'target'
        if (line.find(target) != std::string::npos) {
            // Trova la posizione della stringa 'afterThat' nella riga
            size_t pos = line.find(afterThat);
            if (pos != std::string::npos) {
                // Estrai tutto dopo la stringa 'afterThat'
                std::string valuePart = line.substr(pos + afterThat.size());

                // Converti la parte successiva in un double
                std::istringstream iss(valuePart);
                double chiSquared;
                if (iss >> chiSquared) {
                    file.close();
                    return chiSquared;
                }
            }
        }
    }

    // Se il file non contiene 'Target'
    //std::cerr << "Error: value not found in file " << filePath << std::endl;
    file.close();
    throw runtime_error("value not found in file");
    return 1e30; // non è stato trovato, quindi non è ancora mai stato processato il file, quindi di sicuro il chi quadro calcolato migliora il presente (che non c'è appunto)
}



/*
double input::GetChiSquared(const std::string filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file: " << filePath << std::endl;
        return -1; // codice di errore per file non trovato
    }

    std::string line;
    while (std::getline(file, line)) {
        // Cerca la riga che contiene "chi-squared"
        if (line.find("chi-squared") != std::string::npos) {
            // Trova il primo numero nella riga
            std::istringstream iss(line);
            std::string chiSquaredString;
            double chiSquared;

            // Leggi fino a trovare il valore
            if (std::getline(iss, chiSquaredString, '=') &&
                iss >> chiSquared) {
                file.close();
                return chiSquared;
            }
        }
    }

    // Se il file non contiene "chi-squared"
    std::cerr << "Error: chi-squared not found in file " << filePath << std::endl;
    file.close();
    return 1e30; // non è stato trovato, quindi non è ancora mai stato processato il file, quindi di sicuro il chi quadro calcolato migliora il presente (che non c'è appunto)
}
*/