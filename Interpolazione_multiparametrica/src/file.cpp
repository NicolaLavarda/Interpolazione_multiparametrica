#include "file.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>


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

    // Legge il file riga per riga
    while (std::getline(file, line)) {

        // Sostituisce le virgole con i punti nella riga
        replace(line.begin(), line.end(), ',', '.');

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







/*
void Popolamento_vettori(vector<double>& x, vector<double>& y, vector<double>& sigma_y, string file_dati) {

    ifstream infile(file_dati);

    if (!infile.is_open()) {
        cout << "Unable to open the file" << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(infile, line)) {
        // Sostituisce le virgole con i punti nella riga
        replace(line.begin(), line.end(), ',', '.');

        // Utilizza un istringstream per estrarre i valori
        istringstream iss(line);
        double val_x, val_y, val_sigma_y;

        if (iss >> val_x >> val_y >> val_sigma_y) {
            x.push_back(val_x);
            y.push_back(val_y);
            sigma_y.push_back(val_sigma_y);
        }
    }

    infile.close(); // chiusura del file
}
*/