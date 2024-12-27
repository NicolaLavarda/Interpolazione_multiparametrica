#include "util.h"
#include "input.h"

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <regex>

#include <algorithm>

// Funzione ausiliaria per verificare se una stringa è convertibile in numero
bool isNumber(const std::string& str) {
    // Regex per un numero valido (intero o floating point, anche esponenziale)
    // -> Numeri validi: 
    //      Numero intero (123, -456)
    //      Numero decimale (123.45, -0.67, .89)
    //      Notazione scientifica (1.23e4, -4.56E-7)

    std::regex numberPattern(R"(^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$)");
    return std::regex_match(str, numberPattern);
}


// Funzione ausiliaria per verificare la validità dei dati
bool isDataLine(const std::string& line) {
    std::istringstream iss(line);
    std::string word;
    while (iss >> word) {
        if (!isNumber(word)) {
            return false; // La riga non è valida se contiene una parola che non è un numero
        }
    }
    return true; // La riga è valida solo se tutte le parole sono numeri
}


// Funzione ausiliaria per trovare il minimo (in valore assoluto) di un vettore ignorando zero
double findMinIgnoringZero(const std::vector<double>& vec) {
    // Se il vettore è vuoto, restituisci un valore di default
    if (vec.empty()) return 0;

    // Trova il minimo assoluto
    auto minElementIt = std::min_element(vec.begin(), vec.end(), [](double a, double b) {
        return std::fabs(a) < std::fabs(b);
        });

    // Se il minimo è zero, trova il successivo minimo
    if (std::fabs(*minElementIt) == 0.0) {
        // Cerca il minimo ignorando esattamente zero
        auto secondMinIt = std::min_element(vec.begin(), vec.end(), [](double a, double b) {
            if (std::fabs(a) == 0.0) return false; // Ignora zero
            if (std::fabs(b) == 0.0) return true;  // Priorità a valori diversi da zero
            return std::fabs(a) < std::fabs(b);
            });

        // Se non esistono altri elementi, restituisci zero
        return (secondMinIt != vec.end()) ? *secondMinIt : 0.0;
    }

    // Restituisci il minimo trovato
    return *minElementIt;
}


/*
// NON FUNZIONA E NON CAPISCO PERCHE'
double findMinIgnoringZero(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;

    double minValue = fabs(vec[0]); // Inizializza con un valore
    bool found = false;

    for (const double& val : vec) {
        if (val != 0.0) {
            double absVal = std::fabs(val);
            if (absVal < minValue) {
                minValue = absVal;
                found = true;
            }
        }
    }

    return found ? minValue : 0.0; // Restituisci il minimo, oppure 0.0 se non trovato
}
*/