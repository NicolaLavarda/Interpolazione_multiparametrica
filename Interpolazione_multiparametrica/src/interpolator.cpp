#include "interpolator.h"

#include "exprtk.hpp"
#include "util.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <regex>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <iomanip>


static bool flag_in_exp = false;
static bool flag_in_data = false;

// Costruttore privato
Interpolator::Interpolator() {

    base_order = 10;     // numero più vicino a cui normalizzare i parametri

    // Inizializzo gli ordini di correzione dei parametri a 1 solo la prima volta che creo la prima istanza
    if (order_par.empty())
        order_par.resize(name_par.size(), 1.0);

}


Interpolator* Interpolator::getNewInstance() {
    if (!(flag_in_exp && flag_in_data)) {
        throw std::runtime_error("Accessing new Interpolator istance without primary setting data");
        return nullptr;
    }
    Interpolator& i_generator = Interpolator::getInstance();
    return i_generator.setNewInstance();
}

Interpolator* Interpolator::setNewInstance() {
    Interpolator* instance = new Interpolator();

    instance->setExpression(f_interpolante);      // imposta le variabili uguali a qualsiasi istanza Interpolator
    instance->setData(x, sigma_x, y, sigma_y);    //

    return instance;
}

// Accesso all'istanza permanente (termina alla fine del programma)
Interpolator& Interpolator::getInstance() {
    static Interpolator instance;
    return instance;
}

void Interpolator::setExpression(const std::string& espressione_interpolante) {
    f_interpolante = espressione_interpolante;

    static std::string exp_val = espressione_interpolante;      // La prima volta in assoluto che viene chiamata questa funzione viene definito 'f_interpolante_const' pari a quella effettiva data in input
    f_interpolante_const = exp_val;                 // rimane poi la stessa per qualsiasi istanza a questa classe

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
        throw std::runtime_error("Error in 'setExpression': Interpolating expression registration failed.");
    }

    flag_in_exp = true;
}


void Interpolator::setData(std::vector<double>& x_val, std::vector<double>& sigma_x_val,
                           std::vector<double>& y_val, std::vector<double>& sigma_y_val) {

    x = x_val;
    sigma_x = sigma_x_val;
    y = y_val;
    sigma_y = sigma_y_val;

    x_size = x.size();

    flag_in_data = true;
}


void Interpolator::setOrder(int par_num, double order, std::string& f_interp) {      // ad esempio setOrder(0, 0.001) pone al posto del parametro "a" (numero 0) la stringa "a*0.001"
    
    if (par_num > name_par.size()) throw std::invalid_argument("Cannot setOrder of a parameter out of range of defaults");
    
    std::string parametro = name_par[par_num];

    // Costruisco l'espressione regolare per cercare il parametro come parola intera
    std::regex pattern("\\b" + parametro + "\\b");
    // La sostituzione sarà parametro*order, ad es. "b*1000"
    std::string sostituzione = "(" + parametro + "*" + std::to_string(order) + ")";
    f_interp = std::regex_replace(f_interp, pattern, sostituzione);

    //std::cout << "-> funzione: " << f_interp << std::endl;
    //std::cout << "-> par: " << par_num << "\t ordine: " << std::fixed << std::setprecision(8) << order << std::endl;

}


void Interpolator::normalizeTo10(std::vector<double>& par) {

    /*
    std::cout << std::endl << "Prima:" << std::endl;
    for (int i = 0; i < par.size(); i++)
    {
        std::cout << "par" << i << " " << std::fixed << std::setprecision(8) << par[i] << "\t";
        std::cout << "order" << i << " " << std::fixed << std::setprecision(8) << order_par[i] << std::endl;
    }
    std::cout << std::endl << "Dopo:" << std::endl;
    */

    std::string f_interp = f_interpolante_const;    //prende ogni volta la funzione definita all'inizio dall'utente
    
    int n = par.size();
    for (int i = 0; i < n; ++i) {

        double potenza = std::pow(10.0, -std::round(std::log10(base_order / std::fabs(par[i]))));
        //std::cout << "potenza" << i << " " << std::fixed << std::setprecision(8) << potenza << "\t";
        par[i] /= potenza;
        order_par[i] *= potenza;     // aggiorno il vettore static (visibile uguale tra tutte le istanze) che tiene nota delle potenze assegnate ai parametri
        setOrder(i, order_par[i], f_interp);
    }

    /*
    for (int i = 0; i < par.size(); i++)
    {
        std::cout << "par" << i << " " << std::fixed << std::setprecision(8) << par[i] << "\t";
        std::cout << "order" << i << " " << std::fixed << std::setprecision(8) << order_par[i] << std::endl;
    }
    */
    Interpolator& i_generator_base = Interpolator::getInstance();
    i_generator_base.setExpression(f_interp);
}


void Interpolator::denormalize(std::vector<double>& par) {

    int n = par.size();
    for (int i = 0; i < n; ++i) {
        par[i] *= order_par[i];
        order_par[i] = 1;
    }

    std::string f_interp = f_interpolante_const;    //prende ogni volta la funzione definita all'inizio dall'utente
    Interpolator& i_generator_base = Interpolator::getInstance();
    i_generator_base.setExpression(f_interp);
}


std::vector<double> Interpolator::getParOrder() {
    return order_par;
}


double Interpolator::fChiQuadro(std::vector<double> par) {

    //I parametri basta riassegnarli una volta alla chiamata della funzione, la x invece ovviamente va riassegnata ad ogni iterazione del ciclo for per il calcolo del chi quadro
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    double sum_chi = 0;
    static bool err_x = sigma_x.empty();        // per tutto il programma
    if (err_x)
    {
        for (int i = 0; i < x_size; i++)
        {
            symbol_table.get_variable("x")->ref() = double(x[i]);
            sum_chi += std::pow((y[i] - expression.value()) / sigma_y[i], 2);
        }
    }
    else
    {
        for (int i = 0; i < x_size; i++)
        {
            symbol_table.get_variable("x")->ref() = double(x[i]);
            sum_chi += std::pow((y[i] - expression.value()), 2) / (std::pow(sigma_y[i],2) + std::pow(dfdx(i) * sigma_x[i], 2));
        }
    }

    return sum_chi;
}


double Interpolator::yFunction(std::vector<double> par, double x_i) {
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);

    symbol_table.get_variable("x")->ref() = double(x_i);

    return expression.value();
}


//cerco il valore di x[i] in funzione di y[i] con bisezione anziché calcolare la funzione inversa che mi è impossibile in generale
double Interpolator::xFunction(std::vector<double> par, int i) {
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);
    symbol_table.get_variable("d")->ref() = double(par[3]);
    symbol_table.get_variable("e")->ref() = double(par[4]);



    // Definizione della lambda function (con 'this' cattura correttaente le variabili private della classe)
    auto f = [this, i](double x) -> double {
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
        throw std::invalid_argument("Error in 'xFunction': bisection failed");
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
    throw std::runtime_error("Error in 'xFunction': Maximum number of iterations reached without convergence.");
}


double Interpolator::dfdx(int i) {        // O(h^4)
    static double h = findMinIgnoringZero(x) * 1e-6;
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

    return der;
}


int Interpolator::getDimx() {
    return x_size;
}

// Definizione delle variabili statiche
std::vector<double> Interpolator::order_par;
double Interpolator::base_order;

std::vector<std::string> Interpolator::name_par = { "a","b","c", "d", "e" };

std::string Interpolator::f_interpolante_const;