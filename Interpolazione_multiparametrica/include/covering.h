#ifndef COVERING_H
#define COVERING_H

//#include "interpolator.h"

#include <vector>


class covering {
public:

    covering(std::vector<double> par_best, double chi_quadro_min, double cicle_programms, bool output, bool faster);

    void status(int cicle_val, int k);

    void ricoprimento(std::vector<double>& par, std::vector<double>& par_best, int dimensione, bool is_on_surface);

    void next();

    bool exit(bool ricerca_retta);

    bool end();


private:

    std::vector<double> par_prov;
    std::vector<double> par;

    std::vector<double> par_best;
    double chi_quadro_min;

    double chi_quadro_ciclo_prec;

    std::vector<double> chi_quadro;                         //Vettore chi quadri minimi ("ripulito" ad ogni riesecuzione del programma)
    std::vector<int> livelli = std::vector<int>(100, 0);    //Livelli di ogni ricoprimento (valori diversi da 0 se in quel livello viene migliorato il chi quadro)
    int cicle_programms = 1;                                //numero di cicli di riduzione di 'spostamento'
    int livello = 0;

    std::vector<double> passo;
    double spostamento;

    int cicle;
    int k;

    bool output;
    bool faster;

    // Riferimento all'istanza Singleton
    //Interpolator& i_generator = Interpolator::getInstance();

};

#endif