//#include "exprtk.hpp"
#include "file.h"
#include "interpolation.h"
#include "matrix.h"
#include "interpolating_function.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <functional>

using namespace std;

//Variabili globali
vector<double> x, y, sigma_y;                             //Dati iniziali
double chi_quadro_min = 1e10;                             //Chi quadro minimo assoluto (in ogni istante di tutto il programma)
vector<double> chi_quadro;                                //Vettore chi quadri minimi (ripulito ad ogni riesecuzione del programma)
vector<double> par_best;                                  //Parametri effettivamente stampati poi a schermo
vector<int> livelli(100, 0);                              //Livelli di ogni ricoprimento (valori diversi da 0 se in quel livello viene migliorato il chi quadro)
vector<vector<double>> par_matrix;                        //par_best di ogni ciclo di riesecuzione del programma
int cicle_programms = 0;                                  //Numero di cicli di riesecuzione del programma



double funzione_interpolante(vector<double> x, vector<double> par, int i) {
    //return par[0] * sin(x[i] * par[1]) + par[2];
    return par[0] * x[i] + par[1];
    //return par[0] * exp(x[i] * par[1]);
    //return 1 / sqrt(1 + pow(x[i] / par[0], 2));
}



int main(int argc, char* argv[]) {

    if (argc == 1) {
        cout << "No file name has been entered" << endl;
        exit(EXIT_FAILURE);
    }
    Popolamento_vettori(x, y, sigma_y, argv[1]);

    if (argc == 2) {
        cout << "No parameters value have been entered" << endl;
        exit(EXIT_FAILURE);
    }

    /*
    string more = argv[argc - 1];
    
    if (more == "more" || more == "Y") {
        cout << "Faster? [Y/n] "; cin >> faster;
        //cout << "Confidenza nei parametri (es. 0.1=10%): "; cin >> spostamento;
        cout << "Risultati approssimati? [Y/n] "; cin >> approx;
        if (approx == "Y" || approx == "y") approx_bool = 1;
    }
    */


    //Inizializzazione vettore parametri (ricerca logaritmica automatica o valori attorno cui cercare inseriti in compilazione dall'utente)
    vector<double> par;
    //vector<int> par_string;

    int num_a = 0;
    for (int i = 2; i < argc; i++) {
        try {
            par.push_back(stod(argv[i]));
            //par_string.push_back(0);
        }
        catch (const invalid_argument& e) {
            if (string(argv[i]) != "a")
                break;
            if (string(argv[i]) == "a") {
                par.push_back(0);
                //par_string.push_back(1);
                num_a++;
            }
        }
    }

    /*
    for (int i = 0; i < par.size(); i++)
        cout << par[i] << endl;
    */



    
    //Ricerca automatica logaritmica
    if (num_a < par.size()) {
        vector<double> par_auto = par;
        vector<double> par_def = par;
        ricerca_auto(par, par_auto, par_def, 0);
        par = par_def;
    }

    if (num_a == par.size()) {
        vector<double> par_auto = par;
        ricerca_auto(par, par_auto, 0);
    }

    cout << "Parametri da ricerca automatica: " << endl;
    for (int k = 0; k < par.size(); k++)
    {
        cout << par[k] << "\t";
        //par[k] *= 1.69897;          //log(5)
    }


    double sum_chi = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum_chi += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
    }
    chi_quadro_min = sum_chi;
    cout << endl << "Chi quadro = " << sum_chi << endl;
    
    



    //miglioro un altro po' con algoritmo di bisezione
    vector<double> passo_iniziale;
    for (int i = 0; i < par.size(); i++)
    {
        if (fabs(par[i]) < 5) {
            passo_iniziale.push_back(2);       // (se fosse 2 anziché 2.5) così in bisezione iniziale se par[n]=1 cercherà di migliorare tra -1 e 1 perché vicino a 0 potrebbe minimizzare un po' meglio sul segno sbagliato (no ho cambiato... chi do più fiducia... non lo centor in 0)
            par[i] = (par[i] < 0) ? -1 : 1;
        }
        else
            passo_iniziale.push_back(fabs(par[i] * 5));      // [1.69897=log(5)] così in bisezione iniziale se par[n]=100 cercherà di migliorare tra 50 e 150
        
        //passo_iniziale.push_back(fabs(par[i] * 5));      // [1.69897=log(5)] così in bisezione iniziale se par[n]=100 cercherà di migliorare tra 20 e 500
        //cout << passo_iniziale[i] << endl;
    }
    vector<double> par_iniziali(par.size(), 0);
    algoritmo_bisezione(par, par_iniziali, passo_iniziale, 0);    // Bisezione per migliorare un po' la ricerca automatica logaritmica
    par = par_iniziali;

    double sum_chib = 0;
    for (int i = 0; i < x.size(); i++)
    {
        sum_chib += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
    }

    cout << endl;
    cout << "Parametri automatici iniziali migliorati: " << endl;
    for (int k = 0; k < par.size(); k++)
    {
        cout << par[k] << "\t";
    }
    cout << endl;
    cout << "Chi quadro = " << sum_chib << endl;

    chi_quadro_min = 1e10;      //resetto il chi quadro minimo modificato dalla bisezione per il miglioramento della ricerca automatica logaritmica
    par_best = par;


    //Parametri ulteriori inseriti come ultima cosa in comando di compilazione dall'utente
    bool faster = 0;
    bool approx = 0;
    if (argc > par.size() + 2) {
        string more = argv[par.size() + 2];
        for (char c : more) {
            if (c == 'p')
                approx = 1;
            if (c == 'f')
                faster = 1;
        }
    }

    if (faster) {        //Stampa subito a schermo dei risultati
        risultato1(par_best, approx);
        risultato2(par_best, approx);
    }


    //eseguo la 'main' fino ad una effettiva variazione del chi quadro minore del per mille (0.001)
    double chi_quadro_ciclo_prec = 0;
    cicle_programms = 1;
    double spostamento = 0.1;       // es. s=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni

    while (true) {
        chi_quadro_ciclo_prec = chi_quadro_min;     //da ciclo programma precedente

        vector<double> par_def(par.size(), 0);   //Vettore parametri definitivi stampati a schermo
        vector<double> passo;           // Passi diversi per ogni dimensione:
        int livello = 0;                // Livello iniziale
        if (cicle_programms > 2) spostamento /= cicle_programms * 2;        //aumentando fino ad un '*2' ('spostamento /= cicle_programms * 2') il tempo d'esecuzione è sempre migliorato, ma da valutare bene come ottimizzare il parametro '*2'
        
        for (int i = 0; i < par.size(); i++)
            (fabs(par[i]) > 0.5) ? passo.push_back(fabs(par[i] * spostamento)) : passo.push_back(fabs(2 * spostamento));

        par_matrix.push_back(par_best);     //Faccio arrivare i par_best di ogni intero ciclo di programma alla funzione 'ricoprimento' come utimo vettore del vettore di vettori par_matrix


        //Calcolo parametri
        int cicle = 1;      //primo ciclo di livelli
        for (int k = 0; k < 15; k++) {

            cout << endl << "Livello " << cicle << "." << livello << ":";

            // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
            ricoprimento(par, par_def, passo, livello, 0, false);


            if (faster || livello > 4)       //dal livello 5 in poi basta che un solo livello non migliori il chi quadro per passare al ciclo successivo
                livelli[livello - 1] = 0;

            if ((livello > 1) && (livelli[livello] == 0 && livelli[livello - 1] == 0))      //Per livelli 2, 3 e 4 è necessario che almeno gli ultimi due livelli non migliorino il chi quadro per passare al ciclo successivo
            {
                livello = 0;
                if (cicle > 2)               //Per il secondo ciclo rimane 'spostamento' uguale al primo ciclo
                    spostamento /= 4;
                par = par_best;
                for (int i = 0; i < par.size(); i++)
                    passo[i] = (fabs(par[i]) > 0.5) ? fabs(par[i] * spostamento) : fabs(2 * spostamento);
                cicle++;
            }
            else
                ++livello;      // Espansione del ricoprimento


            //Esco se il chi quadro non migliora più del per mille (0.001)
            if (((cicle > 1) || (livello > 4)) && ((chi_quadro.size() > 1) && (livello > 2)))
            {
                double check = (chi_quadro[chi_quadro.size() - 2] - chi_quadro[chi_quadro.size() - 1]) - chi_quadro[chi_quadro.size() - 1] * 0.001;
                if (check < 0)
                    break;
            }

            //Se è stato migliorato solo una volta e sono già al 6° tentativo tanto vale uscire
            if (chi_quadro.size() == 1 && k > 5)
                break;

            //Se non è migliorato nemmeno una volta il chi quadro e sono al 7° tentativo tanto vale terminare tutto il programma
            if (chi_quadro.size() == 0 && k > 6) {
                cout << endl << "Non e' possibile migliorare ulteriormente il chi quadro" << endl;
                return 0;       //Intanto termino qui (non funziona sta roba di ripartire)
                string restart;
                cout << "Insistere e ripartire? [Y/n] "; cin >> restart;
                if (restart == "Y") {
                    //spostamento /= cicle_programms * 4;
                    spostamento = 0.1;
                    //spostamento = 0.2;         //Riparto ma non da 0.1 (esageratamente largo, dovrei comunque essere arrivato attorno al minimo)
                    cicle_programms = 0;
                    chi_quadro_ciclo_prec = 0;
                    break;
                }
                else
                    return 0;
            }
        }



        //Stampo a schermo i risultati con errori dal chi+1
        risultato1(par_best, approx);


        //Stampo a schermo i risultati con errori da matrice di covarianza
        risultato2(par_best, approx);
        

        //riaggiorno parametri, libero 'chi_quadro' e passo al ciclo di programma successivo
        par = par_best;
        chi_quadro.clear();
        cicle_programms++;


        //Termino l'intero programma se il chi quadro subisce effettivamente un miglioramento minore del per mille (0.001)
        if (cicle_programms > 1 && chi_quadro_ciclo_prec - chi_quadro_min < 0.001)
            break;        
    }
    

    return 0;
}