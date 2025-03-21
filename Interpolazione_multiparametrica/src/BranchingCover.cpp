#include "BranchingCover.h"

#include "interpolator.h"
#include "gradient_descent_algorithm.h"

#include <iostream>
#include <vector>
#include <algorithm>  // per std::max_element
#include <cmath>
#include <stdexcept>

BranchingCover::BranchingCover(std::vector<double> par_best):
    par_size(par_best.size()), level(0), max_level(8), death_rate(1), step(0.2)
{
    // Generazione delle direzioni iniziali (base per le successive)
    int k = par_size;
    div.assign(par_size * 2, std::vector<double> (k, 0.0));
    for (int i = 0; i < par_size * 2; i++) {
        if (i < k) {
            div[i][i] = par_best[i] * step;
        }
        else {
            div[i][i - k] = -par_best[i - k] * step;
        }
    }
}

void BranchingCover::SetModel(const int& max_level_model, const int& death_rate_model, const double& step_model) {
    max_level = max_level_model;
    death_rate = death_rate_model;
    step = step_model;
}

bool BranchingCover::compute(std::vector<double>& par_best) {
    bool improved = false;

    branching(par_best, std::pow(10, max_level));

    // controllo qual è il migliore
    double min_chi = i_generator->fChiQuadro(par_best);
    for (const std::vector<double>& par : par_div)
    {
        double chi_i = i_generator->fChiQuadro(par);
        if (chi_i < min_chi) {
            par_best = par;
            min_chi = chi_i;
            improved = true;
        }
    }
    return improved;
}

void BranchingCover::branching(std::vector<double> par_best, int address) {
    level++;

    // creo 'div' insieme di gradienti che puntano in diverse direzioni
    std::vector<std::vector<double>> div = division_direction(address);

    int n = div.size();
    for (int i = 0; i < n; i++)
    {
        std::vector<double> par_i = par_best;       // 'par_best' punto da cui partire ('par_best' diventa ogni volta il 'par_i' che l'ha generato)
        bool check_kill = go_on(par_i ,div[i]);    // fa un passo nella direzione del parametro i-esimo

        // salvataggio informazioni
        par_div.push_back(par_i);        // avanzato di un passo lo salvo
        div_steps.push_back(div[i]);     // salvo la direzione in cui sono andato
        improvements.push_back(check_kill);
        int ad = register_par_address(i, address);      // associo alla posizione del vettore nella ramificazione l'indice corrispondente in 'par_div' (in questo momento è oviamente l'ultimo essendo appena stato aggiunto)
        
        //std::cout << ad << " -> " << par_i[0] << "\t" << par_i[1] << "\t" << par_i[2] << "  ->  " << i_generator->fChiQuadro(par_i) << std::endl;
        //std::cout << par_i[0] << "\t" << par_i[1] << "\t" << par_i[2] << "\t" << i_generator->fChiQuadro(par_i) << std::endl;
        
        if (level != max_level && !branch_kill(ad))     // se sono alla fine (max_level) oppure non migliora da tempo questo ramo, allora non lo faccio più continuare/dividere
            branching(par_i, ad);   // si divide in più direzioni e riparte
    }
    level--;
}

int BranchingCover::register_par_address(int i, int address) {
    // 'address' che riceve è quello dei 'par' da cui è stato generato
    address += std::pow(10, max_level - level) * (i + 1);       // 12140 significa che è nel 4° ramo del 1° ramo del 2° ramo del 1° ramo (e quindi ad una profondità 4 su 5 che è quella massima)
    
    address_par[address] = par_div.size() - 1;
    //std::cout << "address: " << address << std::endl;
    return address;
}


bool BranchingCover::go_on(std::vector<double>& par_i, std::vector<double> div) {

    double chi_min = i_generator->fChiQuadro(par_i);

    for (int i = 0; i < par_size; i++)
    {
        if (div[i] != 0)
            par_i[i] += div[i];
    }

    return (i_generator->fChiQuadro(par_i) < chi_min) ? true : false;     // 'true' se ha migliorato il chi_quadro
}


bool BranchingCover::branch_kill(int address) {

    if (contains(par_div[address_par[address]]))
        return true;                                // true = devo killare il ramo se per quel punto è già passato un'altro ramo
    
    int check = 0;
    for (int i = 0; i < death_rate; i++)
    {
        int idx = address_par[address];
        if (!improvements[idx])
            check++;    // non è migliorato

        address = mother_address(address);
        if (address == 0)
            break;
    }

    return check == death_rate ? true : false;      // true = devo killare il ramo
}


int BranchingCover::mother_address(int son_address) const {
    int max = std::pow(10, max_level);
    if (son_address == max) return 0; // Caso particolare

    int temp = son_address;
    int pos = 1; // Posizione della cifra da modificare

    while (temp % 10 == 0) {
        temp /= 10;
        pos *= 10; // Sposta la posizione alla cifra successiva
    }

    return son_address - (temp % 10) * pos;
}

int BranchingCover::GetLevel(int address) const {

    int level = max_level; // Posizione della cifra da modificare

    while (address % 10 == 0) {
        address /= 10;
        level--;
    }

    return level;
}


// Funzione per verificare se un vettore è già presente in par_div
bool BranchingCover::contains(const std::vector<double>& par) const {

    int limit = par_div.size() - ignore;        // Considera solo i primi n-'ignore' vettori ('ignore'=2 membro private)

    if (par_div.size() < ignore)
        return false; // Se non ci sono abbastanza vettori in par_div, restituisci false (quindi non killa il ramo in questo caso)

    for (int i = 0; i < limit; i++) {
        if (std::equal(par_div[i].begin(), par_div[i].end(), par.begin(),
                [this](double a, double b) { return std::fabs(a - b) < epsilon; })) {
            return true;
        }
    }
    return false;
}


// ritorna n gradienti differenti per andare in direzioni diverse
std::vector<std::vector<double>> BranchingCover::division_direction(int address) {
    
    std::vector<std::vector<double>> div_n = div;
    
    if (address != std::pow(10, max_level))     // vuol dire che non sono all'inizio (prima divisione delle direzioni)
    {
        std::vector<double> steps_i = div_steps[address_par[address]];  // Vettore passi che ha generato 'par'

        // Trova l'indice dell'elemento non nullo in 'steps_i' (è la direzione da cui sono arrivato)
        int pos = -1;
        for (int i = 0; i < par_size; i++) {
            if (steps_i[i] != 0.0) {
                pos = i;
                break;
            }
        }

        if (pos == -1 && address != std::pow(10, max_level)) {
            throw std::runtime_error("Errore: steps non contiene elementi non nulli.");
        }

        // Generazione delle direzioni (prendo quelle base 'div' e tolgo quella che mi farebbe tornare da dove sono arrivato)
        if (pos != -1)
            div_n.erase(div_n.begin() + (pos + par_size));     // elimino l'elemento della direzione opposta a quella da cui sono arrivato, altrimenti tornerei indietro

        if (improvements[address_par[address]])
            div[pos][pos] *= step * GetLevel(address);      // se prima è migliorato, la direzione da cui sta progredendo ne aumento il passo (probabilmente è la direzione giusta)
    }

    return div_n;
}




BranchingCover::~BranchingCover() {
    delete i_generator;
}