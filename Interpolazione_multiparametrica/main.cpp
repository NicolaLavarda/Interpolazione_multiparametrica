#include "input.h"
#include "automatic_research.h"
#include "ChiSquareMinimizer.h"
#include "file.h"
#include "plot.h"

#include <vector>
#include <string>
#include <map>

int main(int argc, char* argv[]) {

    //------------LETTURA INPUT-------------------------

    std::vector<double> x, sigma_x, y, sigma_y;                    //Dati iniziali
    std::vector<double> par_best;                                  //Parametri effettivamente migliori stampati poi a schermo
    int num_a = 0;                  //numero di parametri automatici

    std::map<std::string, bool> options = {
        {"faster", false},
        {"approx", false},
        {"complex", false},
        {"retta", false},
        {"save", false},
        {"plot", false}
    };

    input(argc, argv, par_best, num_a, options, x, sigma_x, y, sigma_y);
    int x_size = x.size();
    
    //------------RICERCA AUTOMATICA--------------------
    
    if (num_a != 0) {
        parametri_auto(par_best, options["complex"]);
    }

    //------------ESECUZIONE PROGRAMMA------------------
    
    ChiSquareMinimizer Optimizer(options, x_size);
    Optimizer.begin  (par_best);
    Optimizer.compute(par_best);
    Optimizer.end    (par_best);    

    //------------SALVATAGGIO DEI RISULTATI-----------------

    if (options["save"])
    {
        //salvo i risultati nel file di testo
        writeFile(std::string(argv[1]), x, sigma_x, y, sigma_y, par_best, options["approx"], std::string(argv[par_best.size() + 2]));
    }

    //------------PRODUZIONE DEI GRAFICI--------------------

    if (options["plot"])
    {
        PlotGenerator p_generator(par_best, x, sigma_x, y, sigma_y, std::string(argv[1]));
        p_generator.compute_plot_function();
        p_generator.compute_plot_chi_distribution();
        p_generator.compute_plot_Residuals();
        p_generator.save(".png");
    }

    return 0;
}