#ifndef RESULTS_H
#define RESULTS_H

#include "interpolator.h"

#include <vector>
#include <string>

// Chiamare nel programma 'Results(par_best, approx);'
// in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati


// Classe Base
class Results_base {
public:
    Results_base(std::vector<double> par, bool approx_bool, std::ostream& output);

    //Volendo può essere usata da sola con 'Results_base.general_result(val, err, name, approx)'
    void general_result(double valore, double errore, std::string nome, bool arrotondamento);

protected:
    std::vector<double> par;
    double chi_min;
    bool approx_bool;
    std::ostream& out; // Stream di output
    std::vector<char> name_par = { 'a','b','c','d','e' };   // nomi delle variabili

private:
    // Riferimento all'istanza Singleton
    Interpolator& i_generator = Interpolator::getInstance();
};

// Classe intermedia 1, eredita virtualmente da Base
class Result1 : virtual public Results_base {
public:
    Result1();
};

// Classe intermedia 2, eredita virtualmente da Base
class Result2 : virtual public Results_base {
public:
    //Stampo a schermo i risultati con errori da matrice di covarianza
    Result2();
};

// Classe derivata, eredita da entrambe le classi intermedie
class Results : public Result1, public Result2 {
public:
    Results(std::vector<double> par_derived, bool approx_bool, int dim_x, std::ostream& output);
    //chiama in ordine 'Results_base', 'Result1' e 'Result2' ed infine qui dentro al costruttore di 'Results'
    // 
    // Chiamare nel programma 'Results(par_best, approx);'
    // in modo da creare un oggetto temporaneo che restituisca semplicemente a schermo i risultati
};


#endif