#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>


// Funzione ausiliaria per verificare se una stringa è convertibile in numero
bool isNumber(const std::string& str);


// Funzione ausiliaria per verificare la validità dei dati
bool isDataLine(const std::string& line);


// Funzione ausiliaria per trovare il minimo (in valore assoluto) di un vettore ignorando zero
double findMinIgnoringZero(const std::vector<double>& vec);


#endif