#include "linear_mode.h"

#include "covering.h"
#include "bisection_algorithm.h"
#include "interpolator.h"
#include "matrix.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <exception>
#include <algorithm>


using namespace std;


linear_mode::linear_mode(vector<double> par, bool faster, bool complex):
    faster(faster), complex(complex)
{
    vector<double> vec(par.size() - 1, 0);
    m_lin = vec;
    q_lin = vec;

}


//funzione effettiva da usare nella main
void linear_mode::research(vector<double> par_in, bool& errore_lin, bool& ricerca_retta) {

    par_best = par_in;

    //Modalità con interpolazione lineare dei parametri migliorati i primi 3 livelli
    int num_point = 3;  //3 punti per l'interpolazione della retta
    if (faster)
        num_point = 2;  //2 punti per l'interpolazione della retta

    //Riempio il vettore con i parametri migliori di ogni livello per poi interpolarli con una retta n-dimensionale
    if (par_lin.size() <= num_point) {
        if (par_lin.empty())
            par_lin.push_back(par_best);
        else
        {
            if (par_best != par_lin.back())
                par_lin.push_back(par_best);
        }
        /*
        cout << endl << "->";
        for (int i = 0; i < par_best.size(); i++)
            cout << "\t" << "->" << par_best[i];
        cout << endl;
        */
    }

    //Se è stato raggiunto il numero di punti prefissati inizio col metodo della retta
    if (par_lin.size() == num_point) {

        try    //faccio la trasposta di 'par_lin'
        {
            par_lin = trasposta(par_lin);
            //cout << "Matrice parametri da interpolare linearmente:" << endl;
            //stampaMatrice(par_lin);

        }
        catch (const out_of_range& e)
        {
            errore_lin = true;              //Motivo per non usare mai più il metodo della retta
            cerr << "error transposed matrix in 'linear_mode': " << e.what() << endl;
        }
        catch (const invalid_argument& e)
        {
            errore_lin = true;              //Motivo per non usare mai più il metodo della retta
            cerr << "error transposed matrix in 'linear_mode': " << e.what() << endl;
        }


        for (int i = 1; i < par_best.size(); i++)   //fa tutto in funzione di par0 (vedo dopo in funzione di quale è meglio fare, ma da qualche parte devo pur iniziare per capire quali parametri in funzione di quale è meglio procedere)
        {
            try {
                linearFit(par_lin[0], par_lin[i], m_lin[i - 1], q_lin[i - 1]);
            }
            catch (const invalid_argument& e) {
                errore_lin = true;              //Motivo per non usare mai più il metodo della retta
                cerr << "Error in linearFit: " << e.what() << endl;
                return;
            }
            catch (const runtime_error& e) {
                errore_lin = true;              //Motivo per non usare mai più il metodo della retta
                cerr << "Error in linearFit: " << e.what() << endl;
                return;
            }
        }

        //cerco di capire qual è il parametro peggiore e in funzione di quello ricalcolo i q ed m
        int n_worst = 0;
        try
        {
            n_worst = worst_parameter(m_lin, q_lin);       //ora 'm_lin' e 'q_lin' hanno la stessa dimensione di 'par_lin'
        }
        catch (const out_of_range& e)
        {
            errore_lin = true;
            cerr << "metodo retta (worst_parameter) out_of_range: " << e.what() << endl;
            return;
        }
        //cout << endl << "--> " << n_worst << " <--" << endl;

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
        try
        {
            bisezione_lin(par_bisez_lin, m_lin, q_lin, n_worst);

        }
        catch (const runtime_error& e)
        {
            errore_lin = true;              //Motivo per non usare mai più il metodo della retta
            cerr << "errore in 'linear_mode': " << e.what() << endl;
            return;
        }

        double chi_lin_fin = i_generator->fChiQuadro(par_bisez_lin);

        if (chi_lin_fin < i_generator->fChiQuadro(par_best))
        {
            par_best = par_bisez_lin;
            chi_quadro_min = chi_lin_fin;
            ricerca_retta = true;

            if (complex)
            {
                cout << endl << "Parametri bisez_lin migliorati rispetto a par" << n_worst << ":" << endl;
                for (int k = 0; k < par_best.size(); k++)
                    cout << par_bisez_lin[k] << "\t";
                cout << endl << "Chi quadro = " << chi_lin_fin << endl;
            }
        }
        else
        {
            errore_lin = true;
            if (complex)
            {
                cout << "-> non ha funzionato il metodo della retta con par" << n_worst << " (continuo al vecchio modo)" << endl;
            }
        }
    }
}




void linear_mode::linearFit(vector<double>& x, vector<double>& y, double& m, double& q) {
    // Controllo che i vettori abbiano la stessa dimensione e non siano vuoti
    if (x.size() != y.size() || x.empty()) {
        throw invalid_argument("Error in 'linearFit': The vectors x and y must have the same size and not be empty..");
    }

    int n = x.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;

    // Calcolo delle somme necessarie
    for (int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_xx += x[i] * x[i];
    }

    // Calcolo di m e q
    double denominator = n * sum_xx - sum_x * sum_x;
    if (denominator == 0) {
        throw runtime_error("Error in 'linearFit': Null denominator. Interpolation not defined.");
    }

    m = (n * sum_xy - sum_x * sum_y) / denominator;
    q = (sum_y - m * sum_x) / n;
}


void linear_mode::bisezione_lin(vector<double>& par, vector<double> m, vector<double> q, int n_worst) {

    //Prima parte a step abbastanza piccoli e tutti regolari
    double range_min_par_n = par[n_worst] / 2;
    double range_max_par_n = par[n_worst] * 10;

    if (par[n_worst] < 0)
    {
        range_min_par_n = par[n_worst] * 10;
        range_max_par_n = par[n_worst] / 2;
    }

    //cout << range_min_par_n << " < " << par[0] << " < " << range_max_par_n << endl;

    //Blocco ricerca bisezione lungo la retta (ad esempio se par[0]=10 allora va da 5 a 100 in bisezione)
    if (1 == 1) {

        double precisione_bisezione = 0.00001;

        int controllo_par_n = 0;
        vector<double> chi_par_n, par_chi_n;
        double min_chi_par_n = i_generator->fChiQuadro(par) * 2;    //Oh là, ecco così viene perfetto. So per certo che migliore almeno fino a 'i_generator->fChiQuadro(par)' perché è incluso nell'intervallo di ricerca. Gli permetto valori fino al doppio per trovare quello giusto                                       //Non capisco perché non funziona correttamente se metto 'double min_chi_par_n = 1e20;'
        double sec_min_chi_par_n = 1000000;              //Questo non servirebbe nemmeno inizializzarlo
        double posizione_min_chi_par_n = 0;
        double posizione_sec_min_chi_par_n = 1;

        for (double p = range_min_par_n; p < range_max_par_n + fabs(range_max_par_n * 0.1); p += (range_max_par_n - range_min_par_n))   //primi due elementi agli estremi --> uso 'p < range_max_par_n + fabs(range_max_par_n * 0.1)' anzichè 'p <= range_max_par_n' per evitare di controllare valori tra double molto vicini
        {
            par[n_worst] = p;
            // Definisco gli altri parametri in funzione del peggiore ( 'par[n_worst]' ) utilizzando gli 'm' e 'q' secondo par[n] = m_n * par[n_worst] + q_n;
            for (int i = 0; i < par.size(); i++)
            {
                if (i != n_worst)
                    par[i] = m[i] * par[n_worst] + q[i];
            }
            double sum_chi = i_generator->fChiQuadro(par);
            chi_par_n.push_back(sum_chi);
            par_chi_n.push_back(p);
        }
        if (chi_par_n[0] > chi_par_n[1]) {
            posizione_min_chi_par_n = 1;
            posizione_sec_min_chi_par_n = 0;
        }

        while (fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) > precisione_bisezione)
        {

            if (par_chi_n[posizione_min_chi_par_n] > par_chi_n[posizione_sec_min_chi_par_n])
            {
                par[n_worst] = par_chi_n[posizione_sec_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
            }
            else
            {
                par[n_worst] = par_chi_n[posizione_min_chi_par_n] + fabs(par_chi_n[posizione_min_chi_par_n] - par_chi_n[posizione_sec_min_chi_par_n]) / 2;
            }

            // Definisco gli altri parametri in funzione del peggiore ( 'par[n_worst]' ) utilizzando gli 'm' e 'q' secondo par[n] = m_n * par[n_worst] + q_n;
            for (int i = 0; i < par.size(); i++)
            {
                if (i != n_worst)
                    par[i] = m[i] * par[n_worst] + q[i];
            }
            double sum_chi = i_generator->fChiQuadro(par);
            par_chi_n.push_back(par[n_worst]);
            chi_par_n.push_back(sum_chi);

            /*
            cout << endl;
            for (int k = 0; k < par.size(); k++)
                cout << par[k] << "\t";
            cout << " -> ";
            cout << sum_chi << endl;
            */

            //cout << "-> " << min_chi_par_n << endl;
            //cout << "--> " << sum_chi << endl;
            if (sum_chi < min_chi_par_n)
            {
                sec_min_chi_par_n = min_chi_par_n;
                min_chi_par_n = sum_chi;
                //cout << "---> " << min_chi_par_n << endl;
                posizione_sec_min_chi_par_n = posizione_min_chi_par_n;
                posizione_min_chi_par_n = chi_par_n.size() - 1;
            }
            else {
                sec_min_chi_par_n = sum_chi;
                posizione_sec_min_chi_par_n = chi_par_n.size() - 1;
            }

            controllo_par_n++;
            if (controllo_par_n > 100)
                throw runtime_error("Error in 'bisezione_lin': Bisection failed");
        }

    }


    /*
    cout << endl;
    for (int k = 0; k < par.size(); k++)
        cout << par[k] << "\t";
    cout << " -> ";
    cout << i_generator->fChiQuadro(par) << endl;
    */


    // par[n_worst] ora è stato migliorato con "Blocco ricerca bisezione lungo la retta"

    range_min_par_n = par[n_worst] * 0.7;
    range_max_par_n = par[n_worst] * 1.3;

    if (par[n_worst] < 0)
    {
        range_min_par_n = par[n_worst] * 1.3;
        range_max_par_n = par[n_worst] * 0.7;
    }


    double min_chi_lin_passi = 1e10;

    vector<double> par_lin_finali;
    double passo_i = (range_max_par_n - range_min_par_n) / 100;
    for (double p = range_min_par_n; p < range_max_par_n; p += passo_i)
    {
        vector<double> par_lin_passi = par;
        vector<double> passi_lin;
        par_lin_passi[0] = p;
        for (int i = 0; i < par.size(); i++)
        {
            if (i != n_worst) {
                par_lin_passi[i] = m[i] * par_lin_passi[n_worst] + q[i];
                passi_lin.push_back(fabs(par_lin_passi[i]));                    //Il passo (che poi diventa la metà in bisezione) è pari al parametro così va circa par*0.5 < par < par*1.5
            }
            else
                passi_lin.push_back(passo_i * 2);       //per il parametro peggiore lo lascio variare in bisezione di un passo pari a quello appunto del passo lungo la retta
        }
        /*
        cout << "-------" << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par_lin_passi[k] << "\t";
        cout << " -> ";
        cout << i_generator->fChiQuadro(par_lin_passi) << endl;
        */

        //Algoritmo di bisezione
        bisection_algorithm(par_lin_passi, passi_lin, min_chi_lin_passi);

        /*
        cout << endl;
        for (int k = 0; k < par.size(); k++)
            cout << par_lin_passi[k] << "\t";
        cout << " -> ";
        cout << i_generator->fChiQuadro(par_lin_passi) << endl;
        cout << "-------" << endl;
        */

        par_lin_finali = par_lin_passi;
    }
    par = par_lin_finali;

    /*
    cout << endl;
    for (int k = 0; k < par.size(); k++)
        cout << par[k] << "\t";
    cout << " -> ";
    cout << i_generator->fChiQuadro(par) << endl;
    */

}


// Funzione per trovare il parametro "peggiore" -> devo capire qual è il parametro peggiore in funzione del quale migliorare gli altri con 'linear_mode'
int linear_mode::worst_parameter(vector<double>& m_lin, vector<double>& q_lin) {
    int n_worst = 0;

    m_lin.insert(m_lin.begin() + n_worst, 1.0);     // aumento la dimensione aggiungendo all'inizio (visto che all'inizio gli m e q sono in funzione di par0) in modo da avere la dimensione ugale a 'par_lin'
    q_lin.insert(q_lin.begin() + n_worst, 0.0);
    
    //capisco qual è il parametro peggiore ovvero quello con il coefficiente angolare percentuale peggiore
    double prec_criticality = 0;
    double chi_quadro_original = i_generator->fChiQuadro(par_best);    // Valore del chi quadro attuale
    for (int i = 0; i < m_lin.size(); i++)
    {
        auto min_iter = min_element(par_best.begin(), par_best.end());
        double epsilon = fabs(*min_iter * 1e-3); // Piccola perturbazione

        vector<double> perturbed_params = par_best;

        //perturbazione proporzionale alla sensibilità di variazione del parametro rispetto a par0
        perturbed_params[i] += epsilon * m_lin[i];      //lascio che abbia il segno di 'm_lin[i]'

        // Calcola la variazione del chi quadro
        double chi_quadro_perturbed = i_generator->fChiQuadro(perturbed_params);

        //cout << endl << "par" << i << ": \t" << chi_quadro_original << endl;
        //cout << endl << "par" << i << ": \t" << chi_quadro_perturbed << endl;

        // Calcola la derivata numerica rispetto al parametro i-esimo
        double derivative = fabs(chi_quadro_perturbed - chi_quadro_original) / fabs(epsilon);    // ) / fabs(epsilon * m_lin[i])
        //cout << endl << "par" << i << ": \t" << derivative << endl;

        // Calcola la criticità normalizzata
        double criticality = fabs(derivative );         // * par_best[i] / chi_quadro_original
        //cout << endl << "par" << i << ": \t" << criticality << endl;


        if (criticality > prec_criticality)
        {
            prec_criticality = criticality;
            n_worst = i;
        }
    }

    //riaggiorno gli m e q in funzione del parametro peggiore
    for (int i = 0; i < m_lin.size(); i++)
    {
        m_lin[i] /= m_lin[n_worst];
        q_lin[i] -= m_lin[i] * q_lin[n_worst];  //è giusto, so che 'm_lin[i]' si è già riaggiornato
    }

    return n_worst;
}


linear_mode::~linear_mode() {
    delete i_generator;
}
