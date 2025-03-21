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
//#include <functional>
#include <map>

using namespace std;

input::input(int argc, char* argv[])
    : argc(argc), argv(argv) {}


void input::compute(std::map<std::string, bool>& options) {

    if (argc == 1) {
        cerr << "No file name has been entered" << endl;
        help();
        exit(EXIT_FAILURE);
    }
    else if (std::string(argv[1]) == "help")
    {
        help();
        exit(EXIT_FAILURE);
    }

    Interpolator& i_generator = Interpolator::getInstance();

    string filePath = string(argv[1]);
    std::vector<double> x, sigma_x, y, sigma_y;

    if (readFile(filePath, x, sigma_x, y, sigma_y) < 0) {
        help();
        exit(EXIT_FAILURE);
    }
    i_generator.setData(x, sigma_x, y, sigma_y);
    SetData(x, sigma_x, y, sigma_y);

    if (argc == 2) {
        cerr << "No parameters value have been entered" << endl;
        help();
        exit(EXIT_FAILURE);
    }


    //Inizializzazione vettore parametri (ricerca logaritmica automatica o valori attorno cui cercare inseriti in compilazione dall'utente)    
    vector<double> par;
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
                //num_a++;
            }
            else
                break;          //if (string(argv[i]) != "a") break;
        }
    }
    int par_size = par.size();


    if (argc == par_size + 2) {
        cerr << "No interpolating function has been entered" << endl;
        help();
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
    SetInterpolation(par, espressione_interpolante, filePath);

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


void input::SetData(std::vector<double>& x, std::vector<double>& sigma_x,
                    std::vector<double>& y, std::vector<double>& sigma_y) {
    data.x = x;
    data.sigma_x = sigma_x;
    data.y = y;
    data.sigma_y = sigma_y;
}

void input::SetInterpolation(std::vector<double> par, std::string interpolating_function, std::string filePath) {
    interpolation.par = par;
    interpolation.interpolating_function = interpolating_function;
    interpolation.filePath = filePath;
}

input::Data input::GetData() {
    return data;
}

input::Interpolation input::GetInterpolation() {
    return interpolation;
}


bool input::improved(const std::string filePath, std::vector<double> par_best) {
    Interpolator* i_generator = Interpolator::getNewInstance();       // necessario se ho normalizzato i parametri
    return (GetValueFromFile(interpolation.filePath, "chi-squared", "=") > i_generator->fChiQuadro(par_best)) ? true : false;
}


std::vector<double> input::GetParametersFromFile() {
    std::vector<double> value_parameters;
    for (size_t i = 0; i < interpolation.name_parameters.size(); i++)
    {
        try
        {
            value_parameters.emplace_back(GetParameter(interpolation.name_parameters[i]));
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
    for (size_t i = 0; i < interpolation.name_parameters.size(); i++)
    {
        try
        {
            value_errors.emplace_back(GetErrorOfParameter(interpolation.name_parameters[i]));
        }
        catch (const std::runtime_error& e)
        {
            break;
        }
    }
    return value_errors;
}


double input::GetParameter(const std::string parameter) {
    return GetValueFromFile(interpolation.filePath, parameter + " =", "(");
}


double input::GetErrorOfParameter(const std::string parameter) {
    return GetValueFromFile(interpolation.filePath, parameter + " =", "+-");
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




void input::help() {
    std::cout << "\n Usage: ./interpolate <file.txt> <param1> <param2> ... <interpolating_function> [options]\n";
    std::cout << "\nExample:\n";
    std::cout << "  ./interpolate file.txt 1.5 a a \"a*sin(x*b)+c\" plot improve\n";
    std::cout << "\nArguments:\n";
    std::cout << "  <file.txt>               - Input data file containing values to interpolate.\n";
    std::cout << "  <param1> <param2> ...    - Initial values for parameters in the interpolation function.\n";
    std::cout << "                             Use 'a' to indicate parameters that should be determined automatically.\n";
    std::cout << "  <interpolating_function> - Mathematical function to be used for interpolation.\n";
    std::cout << "                             This must be enclosed in double quotes (\"\").\n";
    std::cout << "\nOptional Flags:\n";
    std::cout << "  improve  - Uses previously computed parameters from the input file as starting values.\n";
    std::cout << "  faster   - Speeds up the interpolation process.\n";
    std::cout << "  approx   - Rounds results to the appropriate significant figures.\n";
    std::cout << "  complex  - Displays intermediate steps in the chi-squared minimization.\n";
    std::cout << "  save     - Saves results to the input file only if the chi-squared improves or is not present.\n";
    std::cout << "  save!    - Forces saving of results to the input file regardless of improvement.\n";
    std::cout << "  plot     - Generates a graphical representation of the interpolation.\n";
    std::cout << "\nNotes:\n";
    std::cout << "- Ensure the function syntax follows mathematical conventions and includes parameters defined in input.\n";
    std::cout << "- The names of the parameters used in the interpolating function must be written in alphabetical order (e.g., if using three parameters, they must be 'a', 'b', 'c').\n";
    std::cout << "- Flags should be placed at the end of the command.\n";
    std::cout << "- Use 'help' as the only argument to display this message.\n";
}
