#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "exprtk.hpp"
#include "util.h"
#include <vector>
#include <string>
#include <mutex>

// Interpolator è un Singleton
class Interpolator {
public:
    // Disabilita copia e assegnazione
    Interpolator(const Interpolator&) = delete;
    Interpolator& operator=(const Interpolator&) = delete;

    // Accesso all'unica istanza
    static Interpolator& getInstance();

    // Metodi per impostare i dati
    void setExpression(const std::string& espressione_interpolante);
    void setData(std::vector<double>& x_val, std::vector<double>& sigma_x_val,
                 std::vector<double>& y_val, std::vector<double>& sigma_y_val);


    // Effettiva implementazione della funzione 'fChiQuadro'
    double fChiQuadroImpl(std::vector<double> par);

    // Funzione 'fChiQuadro' thread-safe
    double fChiQuadro(const std::vector<double>& par) {
        std::lock_guard<std::mutex> lock(mutex_); // Protezione del blocco
        return fChiQuadroImpl(par);
    }

    // Funzionalità principali
    double yFunction(std::vector<double> par, double x_i);
    double xFunction(std::vector<double> par, int i);

private:
    // Costruttore privato per prevenire istanziamenti diretti
    Interpolator();

    // Funzione per 'fChiQuadro' quando è necessario considerare anche gli errori in x
    double dfdx(int i);

    // Dati recuperati dal file '.txt'
    std::vector<double> x, sigma_x, y, sigma_y;

    // Membro privato per l'espressione, simboli e parser
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    // Dichiarazione globale per l'espressione, simboli e parser
    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;

    // Variabili per 'setExpression'
    double x_i = 0.5;     //Devo inizializzare in qualche modo le variabili dell'espressione_interpolante (tanto vengono cambiate e riaggiornate ogni volta che viene richiamata la funzione 'f_chi_quadro')
    double a = 1;
    double b = 2;
    double c = 3;
    double d = 4;
    double e = 5;

    std::mutex mutex_; // Mutex per proteggere fChiQuadro

};

#endif

