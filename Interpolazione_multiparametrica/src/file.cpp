#include "file.h"
#include "input.h"      // per accedere alla funzione 'bool isNumber(const string& str)'
#include "results.h"      // per stampare i risultati sul file dati
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath> // per std::fabs
#include <limits> // per std::numeric_limits

// Funzione ausiliaria per verificare la validità dei dati
bool isDataLine(const std::string& line) {          // utilizza la funzione 'bool isNumber(const string& str)' definita in "input.h"
    std::istringstream iss(line);
    std::string word;
    while (iss >> word) {
        if (!isNumber(word)) {
            return false; // La riga non è valida se contiene una parola che non è un numero
        }
    }
    return true; // La riga è valida solo se tutte le parole sono numeri
}

// Funzione ausiliaria per trovare il minimo ignorando zero
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

int readFile(const std::string filePath,
    std::vector<double>& x, std::vector<double>& sigma_x,
    std::vector<double>& y, std::vector<double>& sigma_y) {

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file: " << filePath << std::endl;
        return -1; // codice di errore per file non trovato
    }

    std::string line;
    int numColonne = 0;
    int check = 0;

    // Legge il file riga per riga
    while (std::getline(file, line)) {

        // Sostituisce le virgole con i punti nella riga
        replace(line.begin(), line.end(), ',', '.');

        // Verifica se la riga è vuota oppure è una riga non di dati
        if (line.empty() || !isDataLine(line)) {
            check++;
            continue; // continua la lettura alla riga successiva quando incontra una riga non valida
        }
        if (check>2)    // se trovo più di due righe vuote o senza dati validi mi fermo perché sicuramente non ci sono più dati utili
        {
            return numColonne;
        }

        std::istringstream iss(line);
        std::vector<double> valori;
        double valore;

        // Leggi tutti i valori della riga
        while (iss >> valore)
            valori.push_back(valore);

        if (valori.size() == 3) {
            // Se ci sono 3 colonne: x, y, sigma_y
            x.push_back(valori[0]);
            y.push_back(valori[1]);
            sigma_y.push_back(valori[2]);
        }
        else if (valori.size() == 4) {
            // Se ci sono 4 colonne: x, sigma_x, y, sigma_y
            x.push_back(valori[0]);
            sigma_x.push_back(valori[1]);
            y.push_back(valori[2]);
            sigma_y.push_back(valori[3]);
        }
        else {
            std::cerr << "Error: unsupported number of columns (found "
                << valori.size() << " columns) in file: " << filePath << std::endl;
            return -2; // codice di errore per formato non valido
        }

        // Aggiorna il numero di colonne (costante per tutte le righe)
        if (numColonne == 0) {
            numColonne = valori.size();
        }
        else if (numColonne != valori.size()) {
            std::cerr << "Error: inconsistent number of columns in the file: " << filePath << std::endl;
            return -3; // codice di errore per incoerenza nel file
        }
    }

    file.close();
    return numColonne; // Ritorna il numero di colonne
}



void writeFile(const std::string filePath,
     std::vector<double>& x, std::vector<double>& sigma_x,
     std::vector<double>& y, std::vector<double>& sigma_y,
    std::vector<double> par_best, bool approx, std::string function) {
    
    std::ofstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file: " << filePath << std::endl;
        return;
    }

    // Calcola il minimo di ogni vettore ignorando gli zeri
    auto min_x = findMinIgnoringZero(x);
    auto min_sigma_x = findMinIgnoringZero(sigma_x);
    auto min_y = findMinIgnoringZero(y);
    auto min_sigma_y = findMinIgnoringZero(sigma_y);


    // Determina il formato per ogni vettore
    auto getFormat = [](double minValue) -> std::ios_base::fmtflags {
        return (std::fabs(minValue) < 1e-6) ? std::ios::scientific : std::ios::fixed;
        };

    auto format_x = getFormat(min_x);
    auto format_sigma_x = getFormat(min_sigma_x);
    auto format_y = getFormat(min_y);
    auto format_sigma_y = getFormat(min_sigma_y);


    for (size_t i = 0; i < x.size(); i++) {
        file << std::setprecision((format_x == std::ios::scientific) ? 7 : 10);
        file << (format_x == std::ios::scientific ? std::scientific : std::fixed) << x[i];

        file << "  \t";

        if (!sigma_x.empty()) {
            file << std::setprecision((format_sigma_x == std::ios::scientific) ? 7 : 10);
            file << (format_sigma_x == std::ios::scientific ? std::scientific : std::fixed) << sigma_x[i] << "  \t";
        }

        file << std::setprecision((format_y == std::ios::scientific) ? 7 : 10);
        file << (format_y == std::ios::scientific ? std::scientific : std::fixed) << y[i];

        file << "  \t";

        if (!sigma_y.empty()) {
            file << std::setprecision((format_sigma_y == std::ios::scientific) ? 7 : 10);
            file << (format_sigma_y == std::ios::scientific ? std::scientific : std::fixed) << sigma_y[i];
        }

        file << std::endl;
    }

    /*
    for (int i = 0; i < x.size(); i++)
        file << x[i] << "\t"
             << (sigma_x.empty() ? "" : std::to_string(sigma_x[i]) + "\t")
             << y[i] << "\t"
             << (sigma_y.empty() ? "" : std::to_string(sigma_y[i]) + "\t") << std::endl;
    */


    file << std::endl << std::endl;

    file << "Interpolating function: \t" << function;

    Results(par_best, approx, file);


    file.close();
}