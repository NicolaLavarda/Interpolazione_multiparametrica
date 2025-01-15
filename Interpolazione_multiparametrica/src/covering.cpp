#include "covering.h"

//#include "interpolator.h"
#include "ThreadPool.h"
#include "bisection_algorithm.h"
#include "matrix.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
//#include <string>
#include <algorithm>  // per std::max_element

#include <thread>
#include <chrono>  // Per std::chrono::milliseconds

//using namespace std;


// costruttore di base
covering::covering(std::vector<double> par_best, double chi_quadro_min, double cicle_programms, bool output, bool faster) :
    par_best(par_best), chi_quadro_min(chi_quadro_min), cicle_programms(cicle_programms), output(output), faster(faster)
{

    par_prov = par_best;
    par = par_best;
    chi_quadro_ciclo_prec = chi_quadro_min;

    livello = 0;                         // Livello iniziale

    spostamento = 0.1;       // es. spostamento=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
    if (cicle_programms > 2) spostamento /= cicle_programms * 2;        //aumentando fino ad un '*2' ('spostamento /= cicle_programms * 2') il tempo d'esecuzione è sempre migliorato, ma da valutare bene come ottimizzare il parametro '*2'

    for (int i = 0; i < par.size(); i++)
        passo.push_back(spostamento * ((fabs(par[i]) > 0.5) ? fabs(par[i]) : 2));        //spostamento * 2 se miniore di 0.5 in modo che il rispettivo 'passo' sia pari a 0.1*2=0.2 (che nella funzione diventerà un range semilargo 0.1)

    pool = new ThreadPool();
}

void covering::status(int cicle_val, int k_val) {

    cicle = cicle_val;
    k = k_val;

    if (!output) { std::cout << "\rNumber of program iterations: " << cicle_programms << flush; }
    else          { std::cout << endl << "Livello " << cicle << "." << livello << ":";           }
}

// Funzione per generare i parallelepipedi n-dimensionali sulla "superficie" del livello corrente con passi diversi in ogni dimensione
void covering::ricoprimento(std::vector<double>& par_best, int dimensione, bool is_on_surface) {
    
    

    /*
    if (GetFinish()) {
        std::cout << "enddddddddd" << std::endl;
        return;     // concludo
    }
    */

    if (dimensione == par.size()) {
        if (is_on_surface) {

            //Algoritmo di bisezione --------------------------------------------------------------------
            
            std::vector<double> par_bis = par;                       // 'par' ha le coordinate giuste dell'attuale centro del cubo n-dimensionale in analisi, quindi dev'essere l'argomento della funzione, ma non può essere modificato (cosa che invece farebbe la funzione), quindi gli passo una copia
            std::vector<double> passo_bis = passo;

            auto lambda = [](std::vector<double> par_l, std::vector<double> passo_l) -> std::vector<double> {

                std::vector<double> par_in = par_l;
                std::vector<double> passo_in = passo_l;

                double sum_chi_l = 0;
                // Chiama la funzione originale
                bisection_algorithm(par_in, passo_in, sum_chi_l);

                /*
                for (size_t i = 0; i < par_l.size(); i++)
                    std::cout << par_l[i] << "\t";
                std::cout << sum_chi_l << std::endl;
                */

                std::vector<double> par_p = par_in;
                
                // Crea un nuovo vector con il valore aggiuntivo
                par_p.push_back(sum_chi_l);
                return par_p;
                };

            /*
            if (GetFinish()) {
                std::cout << "enddddddddd" << std::endl;
                return;     // concludo
            }
            */
            static int counter = 0;
            //std::cout << "finish = " << finish << std::endl;
            if (counter < 500)    // se non devo terminare       !GetFinish()
            {
                // Aggiungo al pool di multi-thread come task la chiamata di una lambdafunction che ritorna i parametri migliorati da 'bisection_algorithm(par_bis, passo, sum_chi_p);'
                analysis.emplace_back(pool->enqueue(lambda, par_bis, passo_bis));
                counter++;
            }

            //-------------------------------------------------------------------------------------------
        }
        return;
    }

    // Generiamo tutti gli spostamenti per questa dimensione usando il passo specifico (double)
    for (double delta = -livello * passo[dimensione]; delta <= (livello * passo[dimensione] + 1e-9); delta += passo[dimensione]) {
        int m = 0;
        for (int t = 0; t < par.size(); t++)
        {
            if ((fabs(par[t] - (par_prov[t] - livello * passo[t])) < 1e-9) || (fabs(par[t] - (par_prov[t] + livello * passo[t])) < 1e-9)) m++;
        }
        bool surface_condition = m;
        //bool surface_condition = (fabs(delta - livello * passo[dimensione]) < 1e-9);  // Precisione double circa pari a 1e-9
        par[dimensione] = par_prov[dimensione] + delta;  // Aggiorna la coordinata corrente
        ricoprimento(par_best, dimensione + 1, is_on_surface || surface_condition);
    }
}


void covering::GetResults(std::vector<double>& par_best, const unsigned int time_milliseconds, const unsigned int min_tasksCompleted) {

    std::this_thread::sleep_for(std::chrono::milliseconds(time_milliseconds));

    if (output) std::cout << std::endl << "Ho aspettato " << time_milliseconds << "ms" << std::endl;

    unsigned int tasksCompleted = pool->getTasksCompleted();
    unsigned int tasksCompleted_min = min_tasksCompleted - std::thread::hardware_concurrency();  // Il distruttore di 'pool' fa completare le n=std::thread::hardware_concurrency() tasks in più in ogni caso
    if (tasksCompleted < tasksCompleted_min) {

        unsigned int longer = static_cast<unsigned int>(tasksCompleted_min * time_milliseconds * 1.0 / tasksCompleted - time_milliseconds);
        if (output) std::cout << "It takes " << longer << "ms longer than expected" << std::endl;
        
        std::this_thread::sleep_for(std::chrono::milliseconds(longer));

        tasksCompleted = tasksCompleted_min;
    }


    tasksCompleted = pool->getTasksCompleted() + std::thread::hardware_concurrency();   // prima del 'delete' altrimenti non posso più accedere ovviamente alle funzioni di 'pool'

    delete pool;        // Fa in tempo a fare al massimo n tasks dove n è il numero di threads, quindi considero nel tempo 'longer' che fa in tempo a fare n=std::thread::hardware_concurrency() tasks in più

    analysis.resize(tasksCompleted);    // ridimensione il vettore contenente i 'future' da cui posso recuperare i risultati



    if (output) std::cout << "Results levels analysed: " << std::endl;
    unsigned int ind = 0;
    for (auto& future : analysis)
    {
        auto res = future.get(); // Evita una copia temporanea qui

        if (res.back() < chi_quadro_min)
        {
            chi_quadro_min = res.back();
            par_best.assign(res.begin(), res.end() - 1); // Assegna direttamente evitando un ciclo.
        }

        if (output)
        {
            std::cout << ind + 1 << ".\t"; ind++;
            for (size_t k = 0; k < par_best.size(); k++)
            {
                std::cout << res[k] << "\t";
            }
            std::cout << res.back() << std::endl;
        }
    }

    



    /*
    std::cout << "Risultati..." << std::endl;
    for (auto& future : analysis)
    {
        auto res = future.get(); // Evita una copia temporanea qui.

        std::cout << res.back() << std::endl;

        if (res.back() < chi_quadro_min)
        {
            chi_quadro_min = res.back();
            par_best.assign(res.begin(), res.end() - 1); // Assegna direttamente evitando un ciclo.
        }
    }
    */

    /*
    for (size_t i = 0; i < tasksCompleted; i++)
    {
        std::vector<double> res = analysis[i].get();

        if (res.back() < chi_quadro_min)
        {
            chi_quadro_min = res.back();
            for (size_t k = 0; k < par_best.size(); k++)
                par_best[k] = res[k];
        }

        if (output)
        {
            std::cout << i + 1 << ".\t";
                for (size_t k = 0; k < par_best.size(); k++)
                {
                    std::cout << res[k] << "\t";
                }
            std::cout << res.back() << std::endl;
        }
    }
    */
}


void covering::next() {
    livello++;
}


void covering::SetFinish(bool value) {
    finish.store(value, std::memory_order_relaxed); // Imposta il valore
}

bool covering::GetFinish() const {
    return finish.load(std::memory_order_relaxed); // Legge il valore
}

bool covering::end() {

    /*
    /riaggiorno parametri, libero 'chi_quadro' e 'livelli', e passo al ciclo di programma successivo
    par = par_best;
    livelli.clear();            //livelli viene ripristinato
    livelli.assign(100, 0);     //
    chi_quadro.clear();         //chi_quadro viene ripristinato

    if (cicle_programms == 1)
    {
        //Non posso uscire già al primo ciclo
        return false;
    }

    //Termino l'intero programma se il chi quadro subisce effettivamente un miglioramento minore del per mille (0.001)
    if (cicle_programms > 1 && (chi_quadro_ciclo_prec != chi_quadro_min && chi_quadro_ciclo_prec - chi_quadro_min < 0.001))
        return true;

    // [*] Se non è migliorato nemmeno una volta il chi quadro e ero al 7° tentativo tanto vale terminare tutto il programma
    if (chi_quadro.size() == 0 && k > 6) {
        return true;
    }
    */

    return false;

}

