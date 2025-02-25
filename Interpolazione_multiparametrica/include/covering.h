#ifndef COVERING_H
#define COVERING_H

//#include "interpolator.h"
#include "ThreadPool.h"

#include <vector>
#include <future>
//#include <atomic>


class covering {
public:

    covering(std::vector<double> par_best, double chi_quadro_min, double cicle_programms, bool output, bool faster);

    void status(int cicle_val, int k);

    void ricoprimento(int dimensione, bool is_on_surface);

    void GetResults(std::vector<double>& par_best, const unsigned int time_milliseconds, const unsigned int min_tasksCompleted, const unsigned int time_max);

    void next();


private:

    std::vector<double> par_prov;   //restano fissi per indicare il centro da cui parte il ricoprimento
    std::vector<double> par;        //coordinate che cambiano del centro dell'attuale cubo n-dim in analisi

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

    // opzioni richieste in linea di comando
    bool output;
    bool faster;

    // Pool di multi-thread
    ThreadPool* pool = nullptr;
    std::vector<std::future<std::vector<double>>> analysis;
    unsigned int numThreads = 20;  //std::thread::hardware_concurrency();

    int counter = 0;
    unsigned int max_tasksCompleted = 1000; // valore modificato in 'ChiSquareMinimizer' attraverso 'GetResults'

};

#endif