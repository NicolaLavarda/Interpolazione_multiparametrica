#include "exprtk.hpp"
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
double chi_quadro_min = 10000000;                         //Chi quadro minimo assoluto (in ogni istante di tutto il programma)
vector<double> chi_quadro;                                //Vettore chi quadri minimi (ripulito ad ogni riesecuzione del programma)
vector<double> par_best;                                  //Parametri effettivamente stampati poi a schermo
vector<int> livelli(100, 0);                              //Livelli di ogni ricoprimento (valori diversi da 0 se in quel livello viene migliorato il chi quadro)
vector<vector<double>> par_matrix;                        //par_best di ogni ciclo di riesecuzione del programma
int cicle_programms = 0;                                  //Numero di cicli di riesecuzione del programma


typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

// Dichiarazione globale per l'espressione, simboli e parser
symbol_table_t symbol_table;
expression_t expression;
parser_t parser;

// Variabili globali per i riferimenti a x e par
/*
double* px = nullptr;
double* pa = nullptr;
double* pb = nullptr;
double* pc = nullptr;

double* px = new double;
double* pa = new double;
double* pb = new double;
double* pc = new double;
*/


/*
void setup_expression(const std::string& expression_str) {
    // Registra i riferimenti ai puntatori
    symbol_table.add_variable("x", *px);
    symbol_table.add_variable("a", *pa);
    symbol_table.add_variable("b", *pb);
    symbol_table.add_variable("c", *pc);

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(expression_str, expression)) {
        std::cerr << "Errore nella compilazione dell'espressione: " << expression_str << std::endl;
        throw std::runtime_error("Compilazione fallita.");
    }
}
*/

double x_i = 0;
double a = 1;
double b = 2;
double c = 3;

void setup_expression(const std::string& expression_str) {
    // Registra i riferimenti ai puntatori
    symbol_table.add_variable("x", x_i);
    symbol_table.add_variable("a", a);
    symbol_table.add_variable("b", b);
    symbol_table.add_variable("c", c);

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(expression_str, expression)) {
        std::cerr << "Errore nella compilazione dell'espressione: " << expression_str << std::endl;
        throw std::runtime_error("Compilazione fallita.");
    }
}

double funzione_interpolante(vector<double>& x, vector<double>& par, int i) {
    // Aggiorna i valori di x e dei parametri in symbol_table
    symbol_table.get_variable("x")->ref() = double(x[i]);
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);

    return expression.value();
}


/*
void setup_expression(const std::string& expression_str) {
    // Registra i riferimenti ai puntatori
    symbol_table.add_variable("x", *px);
    symbol_table.add_variable("a", *pa);
    symbol_table.add_variable("b", *pb);
    symbol_table.add_variable("c", *pc);
    cout << *px << "\t" << *pa << "\t" << *pb << "\t" << *pc << endl;
    symbol_table.get_variable("x")->ref() = double(*px);
    symbol_table.get_variable("a")->ref() = double(*pa);
    symbol_table.get_variable("b")->ref() = double(*pb);
    symbol_table.get_variable("c")->ref() = double(*pc);

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(expression_str, expression)) {
        std::cerr << "Errore nella compilazione dell'espressione: " << expression_str << std::endl;
        throw std::runtime_error("Compilazione fallita.");
    }
}

double funzione_interpolante(vector<double>& x, vector<double>& par, int i) {
    // Aggiorna i puntatori con i nuovi valori
    px = &x[i];
    pa = &par[0];
    pb = &par[1];
    pc = &par[2];
    cout << *px << "\t" << *pa << "\t" << *pb << "\t" << *pc << endl;

    // Restituisce il valore calcolato
    return expression.value();
}

*/


/*
// Dichiarazione globale per l'espressione, simboli e parser
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

symbol_table_t symbol_table;
expression_t expression;
parser_t parser;

void setup_expression(string expression_str, vector<double> x, vector<double> par, int i) {
    // Placeholder per le variabili, che verranno aggiornate a ogni chiamata
    symbol_table.add_variable("x", x[i]);
    symbol_table.add_variable("a", par[0]);
    symbol_table.add_variable("b", par[1]);
    symbol_table.add_variable("c", par[2]);

    expression.register_symbol_table(symbol_table);

    if (!parser.compile(expression_str, expression)) {
        cerr << "Errore nella compilazione dell'espressione: " << expression_str << endl;
        throw runtime_error("Compilazione fallita.");
    }
}


double funzione_interpolante(vector<double> x, vector<double> par, int i) {
    // Aggiorna i valori di x e dei parametri in symbol_table
    symbol_table.get_variable("x")->ref() = double(x[i]);
    symbol_table.get_variable("a")->ref() = double(par[0]);
    symbol_table.get_variable("b")->ref() = double(par[1]);
    symbol_table.get_variable("c")->ref() = double(par[2]);

    return expression.value();
}

*/









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

    string more = argv[argc - 1];
    string faster = "";
    bool approx_bool = 1;
    string approx;
    if (more == "more" || more == "Y") {
        cout << "Faster? [Y/n] "; cin >> faster;
        //cout << "Confidenza nei parametri (es. 0.1=10%): "; cin >> spostamento;
        cout << "Risultati approssimati? [Y/n] "; cin >> approx;
        if (approx == "Y" || approx == "y") approx_bool = 1;
    }











    //Inizializzazione vettore parametri (ricerca logaritmica automatica o valori attorno cui cercare inseriti in compilazione dall'utente)
    vector<double> par;
    int par_size = 0;
    if (string(argv[2]) == "a") {
        while ((par_size + 2 < argc) && (string(argv[par_size + 2]) == "a"))
            par_size++;

        vector<double> par_auto(par_size, 0);
        par.resize(par_size, 0);

        //Ricerca automatica logaritmica
        ricerca_auto(par, par_auto, 0);

        double sum_chia = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum_chia += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
        }

        cout << "Parametri automatici iniziali: " << endl;
        for (int k = 0; k < par.size(); k++)
        {
            par[k] *= 10;
            cout << par[k] << "\t";
        }
        cout << endl;
        cout << "Chi quadro = " << sum_chia << endl;

        if (true)
        {
            //miglioro un altro po'
            vector<double> passo_iniziale;
            for (int i = 0; i < par.size(); i++)
            {
                passo_iniziale.push_back(fabs(par[i]*1.5));      // così in bisezione iniziale se par[n]=100 cercherà di migliorare tra 50 e 150
                //cout << passo_iniziale[i] << endl;
            }
            vector<double> par_iniziali(par.size(), 0);
            algoritmo_bisezione(par, par_iniziali, passo_iniziale, 0);    // Bisezione per migliorare un po' la ricerca logaritmica
            
            par = par_iniziali;

            double sum_chib = 0;
            for (int i = 0; i < x.size(); i++)
            {
                sum_chib += pow((y[i] - funzione_interpolante(x, par, i)) / sigma_y[i], 2);
            }

            cout << "Parametri automatici iniziali migliorati: " << endl;
            for (int k = 0; k < par.size(); k++)
            {
                cout << par[k] << "\t";
            }
            cout << endl;
            cout << "Chi quadro = " << sum_chib << endl;
        }
    }
    else
    {
        for (int i = 2; i < argc; i++) {
            try{
                par.push_back(stod(argv[i]));
            }
            catch (const invalid_argument& e){
                break;
            }
        }
    }
    par_best = par;

    if (faster == "Y") {        //Stampa subito a schermo dei risultati
        risultato1(par_best, approx_bool);
        risultato2(par_best, approx_bool);
    }
        


    /*
    *px = 1.0;
    *pa = 2.0;
    *pb = 3.0;
    *pc = 4.0;
    */

    string espressione_interpolante;
    cout << "Inserisci la funzione interpolante (es: a*sin(x*b)+c): ";
    getline(cin, espressione_interpolante);

    // Inizializza l'espressione solo una volta
    setup_expression(espressione_interpolante);
    //setup_expression(espressione_interpolante, x, par, 0);

    /*
    for (int i = 0; i < par.size(); i++)
    {
        cout << par[i] << "\t";
    }
    cout << endl;
    
    for (size_t i = 0; i < x.size(); ++i) {
        //cout << "->" << x[i] << par[0] << par[1] << par[2] << endl;
        double result = funzione_interpolante(x, par, i);
        cout << "Risultato per x[" << i << "] = " << result << endl;
    }
    */

    
    




    //eseguo la 'main' fino ad una effettiva variazione del chi quadro minore del per mille (0.001)
    double chi_quadro_ciclo_prec = 0;
    cicle_programms = 1;
    
    while (true) {
        chi_quadro_ciclo_prec = chi_quadro_min;     //da ciclo programma precedente
        
        vector<double> par_def(par.size(), 0);   //Vettore parametri definitivi stampati a schermo
        vector<double> passo;           // Passi diversi per ogni dimensione:
        int livello = 0;                // Livello iniziale
        double spostamento = 0.1;       // es. s=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
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


            if (faster == "Y" || livello > 4)       //dal livello 5 in poi basta che un solo livello non migliori il chi quadro per passare al ciclo successivo
                livelli[livello - 1] = 0;

            if ((livello > 1) && (livelli[livello] == 0 && livelli[livello - 1] == 0))      //Per livelli 2, 3 e 4 è necessario che almeno gli ultimi due livelli non migliorino il chi quadro per passare al ciclo successivo
            {
                livello = 0;
                if (cicle > 2)
                    spostamento /= 4;           //Per il secondo ciclo rimane 'spostamento' uguale al primo ciclo
                par = par_best;
                for (int i = 0; i < par.size(); i++)
                    passo[i] = (fabs(par[i]) > 0.5) ? fabs(par[i] * spostamento) : fabs(2 * spostamento);
                cicle++;
            }
            else                
                ++livello;      // Espansione del ricoprimento


            //setup_expression(espressione_interpolante, x, par, 0);

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
                return 1;
            }
        }



        //Stampo a schermo i risultati con errori dal chi+1
        risultato1(par_best, approx_bool);


        //Stampo a schermo i risultati con errori da matrice di covarianza
        risultato2(par_best, approx_bool);
        

        //riaggiorno parametri, libero 'chi_quadro' e passo al ciclo di programma successivo
        par = par_best;
        chi_quadro.clear();
        cicle_programms++;


        //Termino l'intero programma se il chi quadro subisce effettivamente un miglioramento minore del per mille (0.001)
        if (cicle_programms > 1 && chi_quadro_ciclo_prec - chi_quadro_min < 0.001)
            break;        
    }
    
    /*
    delete px;
    delete pa;
    delete pb;
    delete pc;
    */

    return 0;
}