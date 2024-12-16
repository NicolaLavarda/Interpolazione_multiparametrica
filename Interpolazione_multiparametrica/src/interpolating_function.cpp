#include "exprtk.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <functional>
#include <algorithm>
#include <stdexcept>

//using namespace std;

extern std::vector<double> x, sigma_x, y, sigma_y;
bool flag = true;

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

// Dichiarazione globale per l'espressione, simboli e parser
symbol_table_t symbol_table;
expression_t expression;
parser_t parser;


double x_i = 0;     //Devo inizializzare in qualche modo le variabili dell'espressione_interpolante (tanto vengono cambiate e riaggiornate ogni volta che viene richiamata la funzione 'f_chi_quadro')
double a = 1;
double b = 2;
double c = 3;
double d = 4;
double e = 5;

void setup_expression(std::string espressione_interpolante) {
    // Registra i riferimenti ai puntatori
    symbol_table.add_variable("x", x_i);
    symbol_table.add_variable("a", a);
    symbol_table.add_variable("b", b);
    symbol_table.add_variable("c", c);
    symbol_table.add_variable("d", d);
    symbol_table.add_variable("e", e);
    symbol_table.add_constants();

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(espressione_interpolante, expression)) {
        throw std::runtime_error("Compilazione fallita.");
    }
}





double dfdx(int i) {        // O(h^4)
    static double h = x.back() * 1e-6;
    //double val_i = expression.value();
    symbol_table.get_variable("x")->ref() = double(x[i] + 2 * h);
    double val_ip2 = expression.value();
    symbol_table.get_variable("x")->ref() = double(x[i] + h);
    double val_ip1 = expression.value();
    symbol_table.get_variable("x")->ref() = double(x[i] - 2 * h);
    double val_im2 = expression.value();
    symbol_table.get_variable("x")->ref() = double(x[i] - h);
    double val_im1 = expression.value();


    symbol_table.get_variable("x")->ref() = double(x[i]);       //riporta "x" al valore precedente
    double der = (-val_ip2 + 8 * val_ip1 - 8 * val_im1 + val_im2) / (12 * h);
    //std::cout << "-----> " << der << std::endl;
    return der;
}



double f_chi_quadro(std::vector<double> par) {
    //I parametri basta riassegnarli una volta alla chiamata della funzione, la x invece ovviamente va riassegnata ad ogni iterazione del ciclo for per il calcolo del chi quadro
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    double sum_chi = 0;

    if (flag)
    {
        for (int i = 0; i < x.size(); i++)
        {
            symbol_table.get_variable("x")->ref() = double(x[i]);
            sum_chi += std::pow((y[i] - expression.value()) / sigma_y[i], 2);
        }
    }
    else
    {
        for (int i = 0; i < x.size(); i++)
        {
            symbol_table.get_variable("x")->ref() = double(x[i]);
            sum_chi += std::pow((y[i] - expression.value()), 2) / (std::pow(sigma_y[i],2) + std::pow(dfdx(i) * sigma_x[i], 2));
        }
    }
    


    return sum_chi;
}


double funzione_interpolante(std::vector<double> par, double x_i) {
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    symbol_table.get_variable("x")->ref() = double(x_i);

    return expression.value();
}




//cerco il valore di x[i] in funzione di y[i] con bisezione anziché calcolare la funzione inversa che mi è impossibile
double x_function(std::vector<double> par, int i) {
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);



    // Definizione della lambda function
    auto f = [i](double x) -> double {
        symbol_table.get_variable("x")->ref() = x;
        return y[i] - expression.value(); // ora 'i' è catturato
    };

    static double min_x = *min_element(x.begin(), x.end());
    static double max_x = *max_element(x.begin(), x.end());


    double a = x[i] - std::fabs(max_x - min_x) * 0.2;       //cerco in un range attorno al valore di x[i] che mi aspetto
    double b = x[i] + std::fabs(max_x - min_x) * 0.2;
    double tollerance = std::fabs(x[i]) * 1e-6;

    //cout << a << "\t" << b << "\t" << f(a) << "\t" << f(b) << endl;

    if (f(a) * f(b) >= 0) {
        throw std::invalid_argument("f(a) e f(b) devono avere segni opposti!");
    }

    double c; // Punto medio
    int max_iter = 100;
    for (int k = 0; k < max_iter; k++) {
        c = (a + b) / 2; // Calcolo del punto medio
        double fc = f(c);

        // Controllo se il risultato è sufficientemente vicino a 0 o se l'intervallo è piccolo
        if (std::fabs(fc) < tollerance || std::fabs(b - a) < tollerance) {
            return c;
        }

        // Aggiornamento dell'intervallo
        if (fc * f(a) < 0) {
            b = c; // Nuovo intervallo: [a, c]
        }
        else {
            a = c; // Nuovo intervallo: [c, b]
        }
    }

    // Se il numero massimo di iterazioni è stato raggiunto
    throw std::runtime_error("Raggiunto il numero massimo di iterazioni senza convergenza.");
}

