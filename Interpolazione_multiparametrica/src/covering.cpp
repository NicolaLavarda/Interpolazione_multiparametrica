#include "covering.h"

//#include "interpolator.h"
#include "ThreadPool.h"
#include "bisection_algorithm.h"
#include "matrix.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>  // per std::max_element
#include <thread>
#include <chrono>  // Per std::chrono::milliseconds


// costruttore di base
covering::covering(std::vector<double> par_best, double chi_quadro_min, double cicle_programms, bool output, bool faster) :
    par_best(par_best), chi_quadro_min(chi_quadro_min), cicle_programms(cicle_programms), output(output), faster(faster)
{

    par_prov = par_best;
    par = par_best;
    //chi_quadro_ciclo_prec = chi_quadro_min;

    livello = 0;                         // Livello iniziale

    spostamento = 0.1;       // es. spostamento=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
    spostamento /= cicle_programms * 2;     // dimezzo lo spostamento ad ogni riesecuzione ('cicle_programms++')

    for (int i = 0; i < par.size(); i++)
        passo.push_back(spostamento * ((fabs(par[i]) > 0.5) ? fabs(par[i]) : 2));        //spostamento * 2 se miniore di 0.5 in modo che il rispettivo 'passo' sia pari a 0.1*2=0.2 (che nella funzione diventerà un range semilargo 0.1)

    pool = new ThreadPool(numThreads);
}


void covering::status(int cicle_val, int k_val) {

    cicle = cicle_val;
    k = k_val;

    if (!output) { std::cout << "\rNumber of program iterations: " << cicle_programms << flush; }
    //else          { std::cout << endl << "Livello " << cicle << "." << livello << ":";           }

}

// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void covering::ricoprimento(int dimensione, bool is_on_surface) {

    if (dimensione == par.size()) {
        if (is_on_surface) {

            //Algoritmo di bisezione --------------------------------------------------------------------
            
            auto lambda = [](std::vector<double> par_l, std::vector<double> passo_l) -> std::vector<double> {

                double sum_chi_l = 0;
                // Chiama l'algoritmo di bisezione
                bisection_algorithm(par_l, passo_l, sum_chi_l);

                /*
                for (size_t i = 0; i < par_l.size(); i++)
                    std::cout << par_l[i] << "\t";
                std::cout << sum_chi_l << std::endl;
                */
                
                // Crea un nuovo vector con il valore aggiuntivo
                par_l.push_back(sum_chi_l);

                return par_l;
            };


            //std::cout << "finish = " << finish << std::endl;
            if (counter < max_tasksCompleted)    // se non devo terminare       !GetFinish()
            {
                // Aggiungo al pool di multi-thread come task la chiamata di una lambdafunction che ritorna i parametri migliorati da 'bisection_algorithm(par_bis, passo, sum_chi_p);'
                analysis.emplace_back(pool->enqueue(lambda, par, passo));
                counter++;
                //std::cout << "\t" << counter << std::endl;
            }

            //-------------------------------------------------------------------------------------------
        }
        return;
    }

    // Generiamo tutti gli spostamenti per questa dimensione usando il passo specifico (double)
    for (double delta = -livello * passo[dimensione]; delta <= (livello * passo[dimensione] + 1e-9); delta += passo[dimensione]) {
        par[dimensione] = par_prov[dimensione] + delta;  // Aggiorna la coordinata corrente
        
        int m = 0;
        for (int t = 0; t < par.size(); t++)
        {
            if ((fabs(par[t] - (par_prov[t] - livello * passo[t])) < 1e-9) || (fabs(par[t] - (par_prov[t] + livello * passo[t])) < 1e-9)) m++;
        }
        bool surface_condition = m;
        
        ricoprimento(dimensione + 1, surface_condition);
    }
}


void covering::GetResults(std::vector<double>& par_best, const unsigned int time_analysis, const unsigned int min_tasksCompleted, const unsigned int time_max) {

    std::this_thread::sleep_for(std::chrono::milliseconds(time_analysis));
    //if (output) std::cout << std::endl << "Ho aspettato " << time_analysis << "ms" << std::endl;
    unsigned int tasksCompleted = pool->getTasksCompleted();
    max_tasksCompleted = static_cast<unsigned int>(time_max * tasksCompleted * 1.0 / time_analysis);    //in funzione del massimo tempo che voglio impiegarci decido di conseguenza quante tasks fare
            //std::cout << "max_tasksCompleted = " << max_tasksCompleted << std::endl;
    unsigned int tasksCompleted_min = ((min_tasksCompleted < max_tasksCompleted) ? min_tasksCompleted : max_tasksCompleted - numThreads) - numThreads;  // Il distruttore di 'pool' fa completare le n=numThreads tasks in più in ogni caso
            //std::cout << "tasksCompleted_min = " << tasksCompleted_min << std::endl;
    // se le tasks fatte sono meno del minimo calcolo il tempo di cui ha ancora bisogno per arrivare al minimo
    if (tasksCompleted < tasksCompleted_min) {

        unsigned int longer = static_cast<unsigned int>(tasksCompleted_min * time_analysis * 1.0 / tasksCompleted - time_analysis);
        if (output) std::cout << "It takes " << longer << "ms longer than expected" << std::endl;
            //std::cout << "Tasks completed: " << pool->getTasksCompleted() << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(longer));

        tasksCompleted = tasksCompleted_min;
    }
    //std::cout << "Tasks completed: " << pool->getTasksCompleted() << std::endl;
    // se le tasks fatte sono più del massimo (tanto meglio) le prendo tutte (il massimo meno 'numThreads' perché n='numThreads' tasks vengono comunque fatte nel 'delete' del pool di threads)
    //tasksCompleted = pool->getTasksCompleted() + numThreads;   // prima del 'delete' altrimenti non posso più accedere ovviamente alle funzioni di 'pool'
    //if (tasksCompleted >= max_tasksCompleted - numThreads) tasksCompleted = max_tasksCompleted - numThreads;

    // Termino il tutto chiamando il 'delete' di 'pool' (facendo 'SetEndAllTasks(false)' faccio in modo che non termini tutte quelle ancora in coda ma solo, come detto, le ultime n='numThreads')
    tasksCompleted = pool->getTasksCompleted();
    pool->SetEndAllTasks(false);      // Faccio in modo che non finisca tutte le tasks in coda ma solo le ultime n, dove n è il numero di threads
    delete pool;        // Fa in tempo a fare al massimo n tasks dove n è il numero di threads, quindi considero nel tempo 'longer' che fa in tempo a fare n=numThreads tasks in più

    analysis.resize(tasksCompleted);    // ridimensione il vettore contenente i 'future' da cui posso recuperare i risultati

    if (output) std::cout << "Tasks completed: " << tasksCompleted << std::endl;

    if (output) std::cout << "Results levels analysed: " << std::endl;
    unsigned int ind = 1;
    for (auto& future : analysis)
    {
        //std::cout << "-> " << ind << std::endl;
        
        auto res = future.get(); // Evita una copia temporanea qui

        if (res.back() < chi_quadro_min)
        {
            chi_quadro_min = res.back();
            par_best.assign(res.begin(), res.end() - 1); // Assegna direttamente evitando un ciclo.

            if (output)
            {
                std::cout << ind << ".\t"; ind++;
                for (size_t k = 0; k < par_best.size(); k++)
                {
                    std::cout << res[k] << "\t";
                }
                std::cout << res.back() << std::endl;
            }
        }
    }

    /*
    std::cout << std::endl;
    for (int i = 0; i < par_best.size(); i++)
    {
        std::cout << "par" << i << " " << par_best[i] << "\t";
    }
    std::cout << std::endl;
    */
}


void covering::next() {
    livello++;
}


