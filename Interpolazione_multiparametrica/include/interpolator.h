#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "exprtk.hpp"
#include "util.h"

#include <vector>
#include <string>
#include <map>
#include <memory>


class Interpolator {
public:

    // Accesso all'istanza permanente (termina alla fine del programma)
    static Interpolator& getInstance();

    // Accesso ad una nuova istanza già settata
    static Interpolator* getNewInstance();

    // Metodi per impostare i dati
    void setExpression(const std::string& espressione_interpolante);
    void setData(std::vector<double>& x_val, std::vector<double>& sigma_x_val,
                 std::vector<double>& y_val, std::vector<double>& sigma_y_val);

    // Metodo per normalizzare l'ordine di grandezza dei parametri
    static void normalizeTo10(std::vector<double>& par);

    // Riparta i parametri al loro ordine di grandezza originale
    static void denormalize(std::vector<double>& par);

    // Metodo static per ottenere gli ordini di grandezza corretti dei parametri
    static std::vector<double> getParOrder();

    // Effettiva implementazione della funzione 'fChiQuadro'
    double fChiQuadro(std::vector<double> par);

    // Funzionalità principali
    double yFunction(std::vector<double> par, double x_i);
    double xFunction(std::vector<double> par, int i);

    // Ottieni il numero di dati
    int getDimx();

    //Interpolator();

private:
    Interpolator();

    // Utilizzata da 'getNewInstance()' per settare le variabili accessibili solo da una precisa istanza (non static)
    //std::unique_ptr<Interpolator> setNewInstance();

    Interpolator* setNewInstance();

    // Funzione per 'fChiQuadro' quando è necessario considerare anche gli errori in x
    double dfdx(int i);

    // Metodo per modificare l'ordine di grandezza di un parametro
    static void setOrder(int par_num, double order, std::string& f_interp);

    // Dati recuperati dal file '.txt'
    std::vector<double> x, sigma_x, y, sigma_y;

    // funzione da interpretare (da non modificare)
    static std::string f_interpolante_const;

    // funzione da interpretare (modificabile per cambiare l'ordine di grandezza delle variabili)
    std::string f_interpolante;

    // Membro privato per l'espressione, simboli e parser
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    // Dichiarazione per l'espressione, simboli e parser
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

    // Numero di dati
    int x_size;

    // Parametri a disposizione
    static std::vector<std::string> name_par;

    // Vettore per tenere traccia degli ordini modificati dei parametri
    static std::vector<double>  order_par;

    // numero più vicino a cui normalizzare i parametri
    static double base_order;

};

#endif

