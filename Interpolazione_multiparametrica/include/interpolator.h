#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "exprtk.hpp"
#include "util.h"
#include <vector>
#include <string>

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

    // Funzionalità principali
    double fChiQuadro(std::vector<double> par);
    double yFunction(std::vector<double> par, double x_i);
    double xFunction(std::vector<double> par, int i);

private:
    // Costruttore privato per prevenire istanziamenti diretti
    Interpolator();

    // Funzione per 'fChiQuadro' quando è necessario considerare anche gli errori in x
    double dfdx(int i);

    // Dati recuperati dal file '.txt'
    std::vector<double> x, sigma_x, y, sigma_y;

};

#endif

