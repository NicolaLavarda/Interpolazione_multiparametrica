#include "exprtk.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

extern vector<double> x, y, sigma_y;
//extern string espressione_interpolante;

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

// Dichiarazione globale per l'espressione, simboli e parser
symbol_table_t symbol_table;
expression_t expression;
parser_t parser;


double x_i = 0;
double a = 1;
double b = 2;
double c = 3;
double d = 4;
double e = 5;

void setup_expression(string espressione_interpolante) {
    // Registra i riferimenti ai puntatori
    symbol_table.add_variable("x", x_i);
    symbol_table.add_variable("a", a);
    symbol_table.add_variable("b", b);
    symbol_table.add_variable("c", c);
    symbol_table.add_variable("d", d);
    symbol_table.add_variable("e", e);

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(espressione_interpolante, expression)) {
        std::cerr << "Errore nella compilazione dell'espressione: " << espressione_interpolante << std::endl;
        throw std::runtime_error("Compilazione fallita.");
    }
}


double funzione_interpolante(vector<double>& x, vector<double>& par, int i) {
    // Aggiorna i valori di x e dei parametri in symbol_table
    symbol_table.get_variable("x")->ref() = double(x[i]);
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    return expression.value();
}



double funzione_interpolante1(vector<double> x, vector<double> par, int i) {
    //return par[0] * sin(x[i] * par[1]) + par[2];
    return par[0] * x[i] + par[1];
    //return par[0] * exp(-x[i] * par[1]);
    //return 1 / sqrt(1 + pow(x[i] / par[0], 2));
}



double f_chi_quadro(vector<double> par) {
    //I parametri basta riassegnarli una volta alla chiamata della funzione, la x invece ovviamente va riassegnata ad ogni iterazione del ciclo for per il calcolo del chi quadro
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    double sum_chi = 0;
    for (int i = 0; i < x.size(); i++)
    {
        symbol_table.get_variable("x")->ref() = double(x[i]);
        sum_chi += pow((y[i] - expression.value()) / sigma_y[i], 2);
    }
    return sum_chi;
}