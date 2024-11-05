#include "exprtk.hpp"       //Per l'interpretazione di una stringa come funzione
#include "file.h"
#include "interpolation.h"
#include "matrix.h"
#include "interpolating_function.h"
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
#include <functional>

using namespace std;

double funzione_interpolante(string expression_str, vector<double> x, vector<double> par, int i) {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x[i]);
    symbol_table.add_variable("a", par[0]);
    symbol_table.add_variable("b", par[1]);
    symbol_table.add_variable("c", par[2]);

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    if (!parser.compile(expression_str, expression)) {
        cerr << "Errore nella compilazione dell'espressione: " << expression_str << endl;
        return 0.0;
    }

    return expression.value();
}

int main() {
    vector<double> x = { 1.0, 2.0, 3.0 };
    vector<double> par = { 1.0, 0.5, 2.0 };

    string espressione_interpolante;
    cout << "Inserisci la funzione interpolante (es: a*sin(x*b)+c): ";
    getline(cin, espressione_interpolante);

    for (size_t i = 0; i < x.size(); ++i) {
        double result = funzione_interpolante(espressione_interpolante, x, par, i);
        cout << "Risultato per x[" << i << "] = " << result << endl;
    }

    return 0;
}
