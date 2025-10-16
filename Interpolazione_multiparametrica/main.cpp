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

    std::map<std::string, bool> options = {
        {"improve", false},     // imposta i parametri di partenza pari a quelli gi� calcolati in precedenza e riportati nel file di input
        {"faster", false},      // fa pi� veloce
        {"approx", false},      // approssima i risultati con le giuste cifre significative
        {"complex", false},     // mostra i "passaggi intermedi" nella minimizzazione del chi quadro
        {"save", false},        // salva i risultati nel file contenente i dati solo se il chi quadro migliora rispetto a quello presente nel file (oppure se non presente)
        {"save!", false},       // salva i risultati nel file contenente i dati in ogni caso
        {"plot", false}         // genera i grafici
    };
    
    
    input Input(argc, argv);
    Input.compute(options);

    input::Data data = Input.GetData();
    input::Interpolation interpolation = Input.GetInterpolation();

    std::vector<double> par_best;

    if (options["improve"])
        par_best = Input.GetParametersFromFile();       // parametri presi dal file
    else
        par_best = interpolation.par;                   // parametri presi da riga di comando

    
    //------------RICERCA AUTOMATICA--------------------
    
    AutomaticResearch Auto(par_best, options["complex"]);
    Auto.beginJob();
    Auto.endJob(par_best);

    //------------ESECUZIONE PROGRAMMA------------------

    ChiSquareMinimizer Optimizer(options);
    Optimizer.begin  (par_best);
    Optimizer.compute(par_best);
    Optimizer.end    (par_best);

    interpolation.par = par_best;

    //------------SALVATAGGIO DEI RISULTATI-----------------

    if (options["save!"] || (options["save"] && Input.improved(interpolation.filePath, par_best)))
    {
        //salvo i risultati nel file di testo
        writeFile(data, interpolation, options["approx"]);
    }

    //------------PRODUZIONE DEI GRAFICI--------------------

    if (options["plot"])
    {
        PlotGenerator p_generator(data, interpolation);
        p_generator.compute_plot_function();
        p_generator.compute_plot_chi_distribution();
        p_generator.compute_plot_Residuals();
        p_generator.save(".png");
    }

    return 0;
}