#include "input.h"
#include "automatic_research.h"
#include "ChiSquareMinimizer.h"
#include "file.h"
#include "plot.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>

int main(int argc, char* argv[]) {

    //------------LETTURA INPUT-------------------------

    std::vector<double> x, sigma_x, y, sigma_y;                    //Dati iniziali
    std::vector<double> par_best;                                  //Parametri effettivamente migliori stampati poi a schermo
    int num_a = 0;                  //numero di parametri automatici

    std::map<std::string, bool> options = {
        {"improve", false},     // imposta i parametri di partenza pari a quelli già calcolati in precedenza e riportati nel file di input
        {"faster", false},      // fa più veloce
        {"approx", false},      // approssima i risultati con le giuste cifre significative
        {"complex", false},     // mostra i "passaggi intermedi" nella minimizzazione del chi quadro
        {"retta", false},       // obbliga ad usare il metodo 'linear_mode' ad ogni 'cicle_programms'
        {"save", false},        // salva i risultati nel file contenente i dati solo se il chi quadro migliora rispetto a quello presente nel file (oppure se non presente)
        {"save!", false},       // salva i risultati nel file contenente i dati in ogni caso
        {"plot", false}         // genera i grafici
    };
    
    input Input(argc, argv);
    Input.compute(par_best, num_a, options, x, sigma_x, y, sigma_y);
    if (options["improve"]) {
        par_best = Input.GetParametersFromFile();
        num_a = 0;  //non eseguo la ricerca automatica
    }
    
    //------------RICERCA AUTOMATICA--------------------
    
    if (num_a != 0) {
        AutomaticResearch Auto(par_best, options["complex"]);
        Auto.beginJob();
        Auto.endJob(par_best);
    }

    //------------ESECUZIONE PROGRAMMA------------------

    ChiSquareMinimizer Optimizer(options);
    Optimizer.begin  (par_best);
    Optimizer.compute(par_best);
    Optimizer.end    (par_best);

    //------------SALVATAGGIO DEI RISULTATI-----------------

    if (options["save!"] || (options["save"] && Input.improved(std::string(argv[1]), par_best)))
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