#include "file.h"
#include "input.h"
#include "automatic_research.h"
#include "covering.h"
#include "linear_mode.h"
#include "matrix.h"
#include "interpolating_function.h"
#include "results.h"
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
#include <exception>
#include <functional>
#include <map>
//#include "TCanvas.h"
//#include "TGraph.h"

using namespace std;

//Variabili globali
vector<double> x, sigma_x, y, sigma_y;                    //Dati iniziali
double chi_quadro_min = 1e30;                             //Chi quadro minimo assoluto (in ogni istante di tutto il programma)
vector<double> chi_quadro;                                //Vettore chi quadri minimi (ripulito ad ogni riesecuzione del programma)
vector<double> par_best;                                  //Parametri effettivamente stampati poi a schermo
vector<int> livelli(100, 0);                              //Livelli di ogni ricoprimento (valori diversi da 0 se in quel livello viene migliorato il chi quadro)
vector<vector<double>> par_matrix;                        //par_best di ogni ciclo di riesecuzione del programma
int cicle_programms = 0;                                  //Numero di cicli di riesecuzione del programma


int main(int argc, char* argv[]) {

    //------------LETTURA INPUT-------------------------

    vector<double> par;
    int num_a = 0;

    std::map<std::string, bool> options = {
        {"faster", false},
        {"approx", false},
        {"complex", false},
        {"retta", false},
        {"plot", false}
    };

    input(argc, argv, par, num_a, options);

    bool faster = options["faster"] ? true : false;
    bool approx = options["approx"] ? true : false;
    bool complex = options["complex"] ? true : false;
    bool retta = options["retta"] ? true : false;
    bool plot = options["plot"] ? true : false;


    //------------RICERCA AUTOMATICA--------------------

    if (num_a != 0)
        parametri_auto(par);
    par_best = par;


    if (faster) {        //Stampa subito a schermo dei risultati
        Results(par_best, approx);
    }

    chi_quadro_min = 1e30;      //resetto il chi quadro minimo modificato dalla bisezione per il miglioramento della ricerca automatica logaritmica, in modo che le bisezioni successive in 'ricoprimento' possano partire anche con valori più alti rispetto all'attuale miglior chi quadro minimo


    //------------ESECUZIONE PROGRAMMA------------------

    //eseguo la "main" (intendo la parte vera e propria del programma) fino ad una effettiva variazione del chi quadro minore del per mille (0.001) -> condizione presente alla fine del ciclo while
    double chi_quadro_ciclo_prec = 0;
    cicle_programms = 1;
    bool errore_lin = false;

    while (true) {  //la condizione di termine è alla fine
        chi_quadro_ciclo_prec = chi_quadro_min;     //da ciclo programma precedente        

        vector<double> par_def(par.size(), 0);   //Vettore parametri definitivi stampati a schermo
        vector<double> passo;                    // Passi diversi per ogni dimensione:
        int livello = 0;                         // Livello iniziale

        double spostamento = 0.1;       // es. spostamento=0.1 allora il primo cubo n-dim di ricoprimento è largo il 10% di ogni parametro nelle rispettive direzioni
        if (cicle_programms > 2) spostamento /= cicle_programms * 2;        //aumentando fino ad un '*2' ('spostamento /= cicle_programms * 2') il tempo d'esecuzione è sempre migliorato, ma da valutare bene come ottimizzare il parametro '*2'
        
        for (int i = 0; i < par.size(); i++)
            passo.push_back(spostamento * ( (fabs(par[i]) > 0.5) ? fabs(par[i]) : 2 ) );        //spostamento * 2 se miniore di 0.5 in modo che il rispettivo 'passo' sia pari a 1 (che nella funzione diventerà un range semilargo 0.5)

        par_matrix.push_back(par_best);     //Faccio arrivare i par_best di ogni intero ciclo di programma alla funzione 'ricoprimento' come utimo vettore del vettore di vettori par_matrix

        //Blocco inizializzazioni per ricerca lungo retta -> metodo della retta (molto più efficente dei "ricoprimenti" quando sono distante dalla soluzione)
        vector<vector<double>> par_lin;              //parametri da interpolare linearmente
        vector<double> m_lin(par.size() - 1, 0);     //coefficienti angolari per parametri da interpolare linearmente
        vector<double> q_lin(par.size() - 1, 0);     //intercette per parametri da interpolare linearmente

        //Calcolo parametri
        int cicle = 1;      //primo ciclo di livelli
        for (int k = 0; k < 15; k++) {

            cout << endl << "Livello " << cicle << "." << livello << ":";
            

            // Generazione dei centri dei parallelepipedi n-dimensionali sulla superficie e stima con bisezione
            ricoprimento(par, par_def, passo, livello, 0, false);   //I migliori parametri sono in 'par_best'            


            if (complex)    // l'analisi è presumibilmente lunga e difficile (forse instabile) quindi voglio vedere i parametri migliorati ad ogni fine di livello [-> volendo, aggiungere all'if '&& cicle_programms == 1']
            {
                cout << "\t -> ";
                for (int k = 0; k < par.size(); k++)
                    cout << par_best[k] << "\t";
                //cout << "--->" << f_chi_quadro(par_best);
            }


            //Modalità con ricerca lungo una retta di miglioramento dei parametri
            bool ricerca_retta = false;     //Cambia in 'true' se il metodo qui sotto della ricerca lungo la retta funziona
            if (par.size() == 1) errore_lin = true;
            if ((!errore_lin || retta) && cicle == 1)       // se da input metto retta=true vuol dire che voglio che venga sempre usato quando possibile il metodo della retta
                linear_mode(par_lin, m_lin, q_lin, errore_lin, ricerca_retta, faster);      //fin tanto che non si hanno almeno un tot di punti prefissati li si raccolgono, poi 'linear_mode' continua effettivamente cercando di migliorare i parametri
            


            // Per passare al cilo successivo e ripartire dal livello 0
            bool fast = false;
            if (faster || livello > 4)       //dal livello 5 in poi basta che un solo livello non migliori il chi quadro per passare al ciclo successivo
                fast = true;        //prima c'era 'livelli[livello - 1] = 0;'

            
            //Passaggio al ciclo successivo, altrimenti semplicemente passaggio al livello successivo
            if ((livello > 1) && (livelli[livello] == 0 && (livelli[livello - 1] == 0 || fast)))      //Per livelli 2, 3 e 4 è necessario che almeno gli ultimi due livelli non migliorino il chi quadro per passare al ciclo successivo
            {
                livello = 0;
                if (cicle > 2)               //Per il secondo ciclo rimane 'spostamento' uguale al primo ciclo
                    spostamento /= 4;
                par = par_best;
                
                livelli.clear();
                livelli.assign(100, 0);     //livelli viene ripristinato

                for (int i = 0; i < par.size(); i++)
                    passo[i] = (fabs(par[i]) > 0.5) ? fabs(par[i] * spostamento) : fabs(2 * spostamento);   //Riduco 'spostamento' per i livelli successivi (all'interno dello stesso 'cicle_programme')
                cicle++;
            }
            else
                ++livello;      // Espansione del ricoprimento ('spostamento' ovviamente rimane uguale, ricopro con cubi n-dimensionali di uguali dimensioni in uno stesso ciclo di ricoprimento)





            // MOTIVI DI USCITA DAL CICLO FOR

            //Per la ricerca lungo la retta
            if (ricerca_retta)
            {
                cout << endl << "esco per la retta" << endl;
                break;
            }

            //Solo se sono dopo il quarto livello (oppure dopo il secondo se non sono più al primo ciclo), allora posso terminare il 'cicle_programme' se il chi_quadro è migliorato più del per mille (in ogni caso questo non basta per terminare l'intero programma)
            if (((cicle > 1) || (livello > 4)) && ((chi_quadro.size() > 1) && (livello > 2)))
            {
                double check = (chi_quadro[chi_quadro.size() - 2] - chi_quadro[chi_quadro.size() - 1]) - chi_quadro[chi_quadro.size() - 1] * 0.001;
                if (check < 0) {
                    cout << endl << "esco chi quadro < 0.001" << endl;
                    break;
                }
            }

            //Se è stato migliorato solo una volta e sono già al 6° tentativo tanto vale uscire
            if (chi_quadro.size() == 1 && k > 5) {
                cout << endl << "esco al sesto tentativo" << endl;
                break;
            }

            //Se non è migliorato nemmeno una volta il chi quadro e sono al 7° tentativo tanto vale terminare tutto il programma
            if (chi_quadro.size() == 0 && k > 6) {
                cout << endl << "Non e' possibile migliorare ulteriormente il chi quadro" << endl;
                goto end_loops;     //torna nella main principale uscendo anche dal while
            }
        }



        //Stampo a schermo i risultati
        Results(par_best, approx);


        //riaggiorno parametri, libero 'chi_quadro' e 'livelli', e passo al ciclo di programma successivo
        par = par_best;
        livelli.clear();            //livelli viene ripristinato
        livelli.assign(100, 0);     //
        chi_quadro.clear();         //chi_quadro viene ripristinato
        cicle_programms++;
        

        //Termino l'intero programma se il chi quadro subisce effettivamente un miglioramento minore del per mille (0.001)
        if (cicle_programms > 1 && (chi_quadro_ciclo_prec != chi_quadro_min && chi_quadro_ciclo_prec - chi_quadro_min < 0.001))
            break;
    }
    end_loops:  //arrivo qui se ho trovato 'goto end_loops;' per uscire dal miglioramento dei parametri non riuscendo a migliorare meglio del per mille il chi quadro


    //------------PRODUZIONE DEI GRAFICI--------------------
    if (!plot)
    {
        string plt;
        cout << endl << "Creare i grafici? [Y/n] "; cin >> plt;      // Sicuro di non volere i grafici?
        plot = (plt == "Y" || plt == "y") ? true : false;
    }
    if (plot)
    {
        cout << endl << "--> Grafici non ancora implementati" << endl;

        /*
        //Creo un TCanvas di Root dove inserisco tutti i grafici necessari
        TCanvas* c1 = new TCanvas("c1", "Canvas dinamico", 1400, 800);

        //Vettore di grafici
        vector<TH2F*> grafici_chi2;
        TMultiGraph* grafico_dati_interpolazione = nullptr;

        //----------------------------------------------------
        vector<double> sigma_par;
        for (int i = 0; i < par.size(); i++)
        {
            sigma_par.push_back(par_best / 1000);       // --> PROVVISORIO <--
        }

        //Creazione dei grafici
        plot(grafici_chi2, grafico_dati_interpolazione, par_best, sigma_par);





        // Decidi cosa fare con il canvas nella main
        c1->Draw();                 // Mostra il canvas
        c1->SaveAs("output.png");   // Salva come PNG
        c1->SaveAs("output.root");  // Salva come file ROOT

        // Pulisci la memoria
        for (auto grafico : grafici) {
            delete grafico;
        }
        delete c1;
        */


        /*
        c1->cd(1);
        grafico_1->Draw("a");
        c1->cd(2);
        //gPad->SetLogx();              //Imposta scala logaritmica sulle x
        grafico_dati_interpolazione->Draw("a");
        legend->Draw();
        c1->cd(1);
        gPad->Pop();

        c1->Update();

        //per i grafici a colori:

                c2->cd(1);
        grafico_1->Draw("colz");
        c2->cd(2);
        grafico_2->Draw("a");
        legend->Draw();
        c2->cd(1);
        gPad->Pop();

        c2->Update();
        */


        



        /*
        TCanvas* c2 = new TCanvas("c2", "c2", 0, 0, 1400, 500);
        c2->Divide(2, 0, 0.01, 0.01);

        TCanvas* c3 = new TCanvas("c3", "c3", 0, 0, 1500, 1000);
        c3->Divide(2, 2, 0.02, 0.01);
        */




        /*
        c3->cd(1);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_1->Draw("colz");
        c3->cd(2);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_2->Draw("colz");
        c3->cd(3);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_3->Draw("colz");
        c3->cd(4);
        gPad->SetLeftMargin(0.11);
        gPad->SetRightMargin(0.12);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        grafico_4->Draw("a");
        legend->Draw();
        c3->cd(1);
        gPad->Pop();


        c3->Update();
        */

    }

    

    return 0;
}