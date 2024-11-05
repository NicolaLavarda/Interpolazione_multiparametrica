#include "file.h"
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

using namespace std;


void Popolamento_vettori(vector<double>& x, vector<double>& y, vector<double>& sigma_y, string file_dati) {

    ifstream infile(file_dati);

    if (!infile.is_open()) {
        cout << "Unable to open the file" << endl;
        exit(EXIT_FAILURE);
    }

    // lettura dei dati dal file
    double val_x, val_y, val_sigma_y;

    while (infile >> val_x >> val_y >> val_sigma_y) {
        x.push_back(val_x);
        y.push_back(val_y);
        sigma_y.push_back(val_sigma_y);
    }
    infile.close();         // chiusura del file
}
