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
//string espressione_interpolante;


int main(int argc, char* argv[]) {

    //--------------------------------------------------
    //------------LETTURA INPUT--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------

    if (argc == 1) {
        cerr << "No file name has been entered" << endl;
        exit(EXIT_FAILURE);
    }
    Popolamento_vettori(x, y, sigma_y, argv[1]);

    if (argc == 2) {
        cerr << "No parameters value have been entered" << endl;
        exit(EXIT_FAILURE);
    }


    //Inizializzazione vettore parametri (ricerca logaritmica automatica o valori attorno cui cercare inseriti in compilazione dall'utente)
    vector<double> par;

    int num_a = 0;
    for (int i = 2; i < argc; i++) {
        try {
            par.push_back(stod(argv[i]));
        }
        catch (const invalid_argument& e) {
            if (string(argv[i]) != "a")
                break;
            if (string(argv[i]) == "a") {
                par.push_back(0);
                num_a++;
            }
        }
    }


    if (argc == par.size() + 2) {
        cerr << "No interpolating function has been entered" << endl;
        exit(EXIT_FAILURE);
    }


    //Parametri ulteriori inseriti come ultima cosa in comando di compilazione dall'utente
    bool faster = false;
    bool approx = false;
    bool complex = false;
    string more = argv[argc - 1];
    for (char c : more) {
        if (c == 'p')
            approx = true;
        if (c == 'f')
            faster = true;
        if (c == 't')
            complex = true;
    }

    //Se esistono parametri validi come ultima cosa in comando di compilazione, allora definisco 'espressione_interpolante' fino al penultimo elemento
    int end = argc;
    if (approx || (faster || complex))    //Almeno uno è presente -> if (approx ||(faster || complex))       if (approx || faster)
        end--;


    //Inizializzazione e interpretazione della funzione interpolante (sommo le stringhe in caso ci fossero degli spazi in mezzo alla funzione digitata dall'utente)
    string espressione_interpolante = "";
    for (int i = par.size() + 2; i < end; i++)
    {
        espressione_interpolante += argv[i];
    }
    // Inizializza l'espressione solo una volta
    setup_expression(espressione_interpolante);



    //--------------------------------------------------
    //------------RICERCA AUTOMATICA--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------

    if (num_a != 0)
    {
        //Ricerca automatica logaritmica
        vector<double> par_auto = par;
        vector<double> par_def = par;
        ricerca_auto(par, par_auto, par_def, 0);

        cout << "Parametri da ricerca automatica:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        chi_quadro_min = f_chi_quadro(par);
        cout << endl << "Chi quadro = " << chi_quadro_min << endl;


        //Cerco di capire se è meglio 'par' o 'par_auto' con algoritmo di bisezione
        vector<double> passo_i1;
        vector<double> passo_i2;

        vector<double> par_i1(par.size(), 0);
        vector<double> par_i2(par.size(), 0);

        algoritmo_bisezione(par,      par_i1, passo_i1, 0);
        algoritmo_bisezione(par_auto, par_i2, passo_i1, 0);

        par = (f_chi_quadro(par_i1) < f_chi_quadro(par_i2)) ? par_i1 : par_i2;
        par_best = par;

        cout << "Parametri automatici iniziali migliorati:" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par[k] << "\t";

        chi_quadro_min = f_chi_quadro(par);
        cout << endl << "Chi quadro = " << chi_quadro_min << endl;


    }
    par_best = par;


    if (faster) {        //Stampa subito a schermo dei risultati
        risultato1(par_best, approx);
        risultato2(par_best, approx);
    }

    chi_quadro_min = 1e10;      //resetto il chi quadro minimo modificato dalla bisezione per il miglioramento della ricerca automatica logaritmica



    //--------------------------------------------------
    //------------ESECUZIONE PROGRAMMA--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------


    //eseguo la 'main' fino ad una effettiva variazione del chi quadro minore del per mille (0.001) -> condizione presente alla fine del ciclo while
    double chi_quadro_ciclo_prec = 0;
    cicle_programms = 1;
    double spostamento = 0.1;       // es. spostamento=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
    bool errore_lin = false;

    while (true) {  //la condizione di termine è alla fine
        chi_quadro_ciclo_prec = chi_quadro_min;     //da ciclo programma precedente

        vector<double> par_def(par.size(), 0);   //Vettore parametri definitivi stampati a schermo
        vector<double> passo;                    // Passi diversi per ogni dimensione:
        int livello = 0;                         // Livello iniziale
        if (cicle_programms > 2) spostamento /= cicle_programms * 2;        //aumentando fino ad un '*2' ('spostamento /= cicle_programms * 2') il tempo d'esecuzione è sempre migliorato, ma da valutare bene come ottimizzare il parametro '*2'
        
        for (int i = 0; i < par.size(); i++)
            passo.push_back(spostamento * ( (fabs(par[i]) > 0.5) ? fabs(par[i]) : 2 ) );        //spostamento * 2 se miniore di 0.5 in modo che il rispettivo 'passo' sia pari a 1 (che nella funzione diventerà un range semilargo 0.5)

        par_matrix.push_back(par_best);     //Faccio arrivare i par_best di ogni intero ciclo di programma alla funzione 'ricoprimento' come utimo vettore del vettore di vettori par_matrix

        //Blocco inizializzazioni per ricerca lungo retta
        vector<vector<double>> par_lin;              //parametri da interpolare linearmente
        vector<double> m_lin(par.size() - 1, 0);     //coefficienti angolari per parametri da interpolare linearmente
        vector<double> q_lin(par.size() - 1, 0);     //intercette per parametri da interpolare linearmente

        //Calcolo parametri
        int cicle = 1;      //primo ciclo di livelli
        for (int k = 0; k < 15; k++) {

            cout << endl << "Livello " << cicle << "." << livello << ":";

            // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
            ricoprimento(par, par_def, passo, livello, 0, false);   //I migliori parametri sono in 'par_best'
            chi_quadro.push_back(f_chi_quadro(par_best));   //Alla fine del ricoprimento riempio con l'ultimo chi quadro migliore


            if (complex)    // l'analisi è presumibilmente lunga e difficile (forse instabile) quindi voglio vedere i parametri migliorati ad ogni fine di livello [-> volendo, aggiungere all'if '&& cicle_programms == 1']
            {
                cout << "\t -> ";
                for (int k = 0; k < par.size(); k++)
                    cout << par_best[k] << "\t";
                //cout << "--->" << f_chi_quadro(par_best);
            }


            //Solo nella prima esecuzione del programma
            bool ricerca_retta = false;     //Cambia in 'true' se il metodo qui sotto della ricerca lungo la retta funziona
            if (!errore_lin)           // se non ci sono mai stati errori nella ricerca lungo la retta, oppure l'ultima volta non ha funzionato, allora non ci provo più, anche perché richiede un costo computazionale non indifferente e probabilmente se non ha funzionato una volta non lo farà nemmeno quelle dopo ( prima usavo 'cicle_programms <= 4')
            {
                //Modalità con interpolazione lineare dei parametri migliorati i primi 3 livelli
                int level_lin = 3;  //4 punti per l'interpolazione della retta
                if (faster)
                    level_lin = 2;  //3 punti per l'interpolazione della retta

                if (cicle == 1 && livello <= level_lin)
                    par_lin.push_back(par_best);

                bool err2 = false;
                if (par_lin[0]== par_lin[1])
                {
                    errore_lin = true;
                    err2 = true;
                }

                if ((!err2 && par.size() > 1) && (cicle == 1 && livello == level_lin)) {
                    par_lin = trasposta(par_lin);
                    /*
                    cout << "Matrice parametri da interpolare linearmente:" << endl;
                    stampaMatrice(par_lin);
                    */
                    for (int i = 1; i < par.size(); i++)
                    {
                        try {
                            linearFit(par_lin[0], par_lin[i], m_lin[i - 1], q_lin[i - 1]);
                        }
                        catch (const exception& e) {
                            errore_lin = true;
                            cerr << "Errore in linearFit: " << e.what() << endl;
                            break;
                        }
                    }

                    /*
                    cout << endl << "gli m sono: " << endl;
                    for (int k = 0; k < m_lin.size(); k++)
                        cout << m_lin[k] << "\t";
                    cout << endl;

                    cout << endl << "i q sono: " << endl;
                    for (int k = 0; k < q_lin.size(); k++)
                        cout << q_lin[k] << "\t";
                    cout << endl;
                    */

                    vector<double> par_bisez_lin = par_best;
                    bisezione_lin(par_bisez_lin, m_lin, q_lin);

                    cout << endl << "Parametri bisez_lin:" << endl;
                    for (int k = 0; k < par.size(); k++)
                        cout << par_bisez_lin[k] << "\t";

                    double chi_lin_fin = f_chi_quadro(par_bisez_lin);
                    cout << endl << "Chi quadro = " << chi_lin_fin << endl;

                    if (chi_lin_fin < f_chi_quadro(par_best))
                    {
                        par_best = par_bisez_lin;
                        chi_quadro_min = chi_lin_fin;
                        ricerca_retta = true;
                    }
                    else
                    {
                        errore_lin = true;
                        cout << "-> non ha funzionato il metodo della retta (continuo al vecchio modo)" << endl;
                    }


                }
            }
            


            // Per passare al cilo successivo e ripartire dal livello 0
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






            // MOTIVI DI USCITA DAL CICLO FOR

            //Per la ricerca lungo la retta
            if (ricerca_retta)
            {
                ricerca_retta = false;
                cout << "esco per la retta" << endl;
                break;
            }

            //Esco se il chi quadro non migliora più del per mille (0.001)
            if (((cicle > 1) || (livello > 4)) && ((chi_quadro.size() > 1) && (livello > 2)))
            {
                double check = (chi_quadro[chi_quadro.size() - 2] - chi_quadro[chi_quadro.size() - 1]) - chi_quadro[chi_quadro.size() - 1] * 0.001;
                if (check < 0) {
                    cout << "esco chi quadro < 0.001" << endl;
                    break;
                }
            }

            //Se è stato migliorato solo una volta e sono già al 6° tentativo tanto vale uscire
            if (chi_quadro.size() == 1 && k > 5) {
                cout << "esco al sesto tentativo" << endl;
                break;
            }

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