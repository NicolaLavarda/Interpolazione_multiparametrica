#include "plot.h"
#include "interpolating_function.h"
#include "chi_square.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <set>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h> // se dovessi usare funzioni matematiche di ROOT
#include <TLine.h>      // Per le linee
#include <TColor.h>     // Per i colori predefiniti
#include <TString.h>    // Per usare Form()


using namespace std;


PlotGenerator::PlotGenerator(vector<double> par, vector<double> x, vector<double> sigma_x, vector<double> y, vector<double> sigma_y, string file_name):
    par(par), x(x), y(y), sigma_x(sigma_x), sigma_y(sigma_y)
{
    base_name = file_name.erase(file_name.size() - 4, 4) + "_";      // rimuove '.txt' dal nome del file dati e lo usa per nominare i file dei grafici
}


void PlotGenerator::compute_plot_function() {

    //Creazione del grafico con i dati e la funzione interpolante
    plot_function(grafico_dati_interpolazione, par);

    //primo grafico corrispondente al plot dei dati con funzione interpolante
    Canvas_grafici.push_back(new TCanvas("c1", "c1", 800, 600));
    Canvas_grafici[0]->cd();
    grafico_dati_interpolazione->Draw("a");
    Canvas_grafici[0]->Update();
}

void PlotGenerator::compute_plot_chi_distribution() {

    //Calcolo gli errori sui parameteri per capire quanto grandi plottare i grafici sulle distrubuzioni del chi quadro
    vector<double> sigma_par;
    double chi_piu_uno_val = 1; if (false) chi_piu_uno_val = chi_quadro_piu_uno(par.size());    //metto if(true) se non voglio l'errore sul singolo parametro ma calcolato correttamente per più parametri con la funzione 'chi_quadro_piu_uno'
    double range = 0.20;   //cerco l'errore al chi+1 entro +-20% del valore del parametro (in caso fosse più grande la funzione allarga la ricerca)
    for (int i = 0; i < par.size(); i++)
        sigma_par.push_back(chi_piu_uno_par_n(par, i, chi_piu_uno_val, range));

    //Creazione dei grafici di distribuzione del chi quadro
    if (par.size() == 2 || par.size() == 3)
        plot_chi_distribution(grafici_chi2, par, sigma_par);

    //grafici distribuzioni del chi_quadro
    if (par.size() == 2 || par.size() == 3)     //solo se ci sono 2 o 3 parametri (altrimenti non ha senso fare questo tipo di grafico)
    {
        for (int i = 0; i < grafici_chi2.size(); i++)
        {
            string name_canvas = "c" + to_string(i + 2);
            Canvas_grafici.push_back(new TCanvas(name_canvas.c_str(), name_canvas.c_str(), 800, 600));

            canvas_chi_distribution(Canvas_grafici[i + 1], grafici_chi2[i], par, sigma_par, i);    // aggiungo le linee a 1 sigma e salvo i grafici
        }
    }
}

void PlotGenerator::compute_plot_Residuals() {

    //Creazione dei grafici dei residui (residui x, residui y e residui xy)
    plot_residui(grafici_residui, par);

    //grafici dei residui
    if (sigma_y.size() != 0)
    {
        for (int i = 0; i < grafici_residui.size(); i++)
        {
            string name_canvas = "c" + to_string(i + 2 + grafici_chi2.size());
            Canvas_grafici.push_back(new TCanvas(name_canvas.c_str(), name_canvas.c_str(), 800, 600));

            canvas_residui(Canvas_grafici[i + 1 + grafici_chi2.size()], grafici_residui[i], i);
        }
    }
}

void PlotGenerator::save(string format_file) {

    if (isValidFormat(format_file))
        format = format_file;           //Altrimenti 'format' è di default ".png"

    for (int i = 0; i < Canvas_grafici.size(); i++)
    {
        string file_name_canvas = base_name + "c" + to_string(i + 1) + format;
        Canvas_grafici[i]->SaveAs(file_name_canvas.c_str());
    }

}


//Creazione griglia valori chi_quadro al variare dei parametri (mappe colore)
vector<vector<double>> PlotGenerator::calcolo_matrice_chi_quadri(vector<double> par, vector<double> sigma_par,double num_sigma, double punti_per_parametro, int primo_par, int secondo_par) {

    vector<double> par_i, par_f, passo_par;

    for (int i = 0; i < par.size(); i++)
    {
        par_i.push_back(par[i] - num_sigma * sigma_par[i]);
        par_f.push_back(par[i] + num_sigma * sigma_par[i]);
        passo_par.push_back((par_f[i] - par_i[i]) / punti_per_parametro);
        //cout << par_i[i] << "\t" << par_f[i] << "\t" << passo_par[i] << endl;
    }

    //Creazione griglia scansione completa definitiva da graficare --------------------------------------------------------------------------

    vector<vector<double>> p1p2chi2(3);     // è un vettore formato da 3 vettori dove il primo indica il primo_par, il secondo il secondo_par e il terzo i relativi chi_quadri

    for (double i = par_i[primo_par]; i <= par_f[primo_par] + fabs(par_f[primo_par]) * 0.00000001; i += passo_par[primo_par])
    {
        for (double j = par_i[secondo_par]; j <= par_f[secondo_par] + fabs(par_f[secondo_par]) * 0.00000001; j += passo_par[secondo_par])
        {
            p1p2chi2[0].push_back(i);
            p1p2chi2[1].push_back(j);
            vector<double> par_chi = par;
            par_chi[primo_par] = i;
            par_chi[secondo_par] = j;
            p1p2chi2[2].push_back(f_chi_quadro(par_chi));
            //cout << i << "\t" << j << "\t" << f_chi_quadro(par_chi) << endl;
        }
    }

    return p1p2chi2;
}


//Grafico dati e funzione interpolante
void PlotGenerator::plot_function(TMultiGraph*& grafico_dati_interpolazione, vector<double> par) {

    gStyle->SetOptStat(0);
    TGraphErrors* dati_plot = new TGraphErrors();
    int dati_size = x.size();
    for (int i = 0; i < dati_size; i++) {
        double err_x = (i < sigma_x.size()) ? sigma_x[i] : 0;
        double err_y = (i < sigma_y.size()) ? sigma_y[i] : 0;

        dati_plot->SetPoint(i, x[i], y[i]);
        dati_plot->SetPointError(i, err_x, err_y);
    }

    dati_plot->SetTitle("Grafico interpolazione");
    dati_plot->GetXaxis()->SetTitle("asse x");
    dati_plot->GetYaxis()->SetTitle("asse y");
    dati_plot->SetMarkerStyle(8);
    dati_plot->SetMarkerSize(0.7);
    dati_plot->SetLineWidth(1);


    TGraph* funzione_interpolante_dati = new TGraph();

    int num_point = 100;     //numero di punti con cui approssimare graficamente la funzione interpolante
    vector<double> vecx_funz_interpolante;

    // minimo e massimo del vettore x
    auto min_iter = min_element(x.begin(), x.end());
    auto max_iter = max_element(x.begin(), x.end());
    double min_val = *min_iter;
    double max_val = *max_iter;
    double step = (max_val - min_val) / num_point;

    for (double i = min_val - step; i < max_val + step; i += step) {
        vecx_funz_interpolante.push_back(i);
    }

    for (int i = 0; i < vecx_funz_interpolante.size(); i++)
    {
        //cout << vecx_funz_interpolante[i] << "\t" << funzione_interpolante(par, vecx_funz_interpolante[i]) << endl;
        funzione_interpolante_dati->SetPoint(i, vecx_funz_interpolante[i], funzione_interpolante(par, vecx_funz_interpolante[i]));        // 'funzione_interpolante(std::vector<double> par, double x)' funzione definita in #include "interpolating_function.h"
    }

    funzione_interpolante_dati->SetLineColor(kRed);
    funzione_interpolante_dati->SetLineWidth(1);

    auto legend = new TLegend(0.7, 0.85, 0.95, 0.95);
    legend->AddEntry(dati_plot, "Dati sperimentali", "p");
    legend->AddEntry(funzione_interpolante_dati, "fit globale", "l");


    grafico_dati_interpolazione = new TMultiGraph("Grafico interpolazione", "Grafico interpolazione");
    grafico_dati_interpolazione->Add(funzione_interpolante_dati, "c");    //mette la funzione interpolante con i dati trovati
    grafico_dati_interpolazione->Add(dati_plot, "p");     //mette i punti dei dati con le loro barre d'errore

}

//Grafici distribuzioni chi quadro
void PlotGenerator::plot_chi_distribution(vector<TH2F*>& grafici_chi2, vector<double> par, vector<double> sigma_par) {

    if (par.size() < 4)     //solo se ci sono 2 o 3 parametri (altrimenti non ha senso fare questo tipo di grafico)
    {
        //nome dei grafici
        vector<string> nome_grafico;
        nome_grafico.push_back("Parametri_a-b");
        nome_grafico.push_back("Parametri_b-c");
        nome_grafico.push_back("Parametri_c-a");

        for (int g = 0; g < 3; g++)
        {
            //cout << "grafico " << g + 1 << endl;
            double num_sigma = 3;
            double punti_per_parametro = 40;
            int primo_par = g;                              // scorro le combianzioni 01 - 12 - 20
            int secondo_par = (g + 1 == 3) ? 0 : g + 1;     //
            vector<vector<double>> p1p2chi2 = calcolo_matrice_chi_quadri(par, sigma_par, num_sigma, punti_per_parametro, primo_par, secondo_par);

            double par_1i = par[primo_par] - num_sigma * sigma_par[primo_par];
            double par_1f = par[primo_par] + num_sigma * sigma_par[primo_par] + fabs(par[primo_par]) * 0.00000001;
            double par_2i = par[secondo_par] - num_sigma * sigma_par[secondo_par];
            double par_2f = par[secondo_par] + num_sigma * sigma_par[secondo_par] + fabs(par[secondo_par]) * 0.00000001;

            auto min_chi = *min_element(p1p2chi2[2].begin(), p1p2chi2[2].end());
            auto max_chi = *max_element(p1p2chi2[2].begin(), p1p2chi2[2].end());

            auto grafico_1 = new TH2F(
                nome_grafico[g].c_str(),              // Nome dell'istogramma
                nome_grafico[g].c_str(),              // Titolo dell'istogramma
                punti_per_parametro,          // Numero di bin sull'asse x
                par_1i,                       // Estensione minima sull'asse x
                par_1f,                       // Estensione massima sull'asse x
                punti_per_parametro,          // Numero di bin sull'asse y
                par_2i,                       // Estensione minima sull'asse y
                par_2f                        // Estensione massima sull'asse y
            );

            grafico_1->SetStats(0);	//Rimuovo il box impostato di defaut da Root

            // Riempi l'istogramma con i parametri e i chi_quadri associati
            for (int i = 0; i < p1p2chi2[0].size(); i++) {
                int bin = grafico_1->FindBin(p1p2chi2[0][i], p1p2chi2[1][i]); // Trova il bin corrispondente alle coordinate (x, y)
                grafico_1->SetBinContent(bin, p1p2chi2[2][i]); // Imposta il valore associato al bin (valore del chi_quadro)
            }

            grafico_1->SetContour(1000);        //quantifica quanto viene "sfumata" la colorazione del grafico

            grafico_1->GetXaxis()->SetTitle("Asse x");
            grafico_1->GetYaxis()->SetTitle("Asse y");
            grafico_1->GetZaxis()->SetTitle("Chi Quadri");

            grafico_1->SetMinimum(min_chi * 0.999999);
            grafico_1->SetMaximum(max_chi * 1.000001);

            grafici_chi2.push_back(grafico_1);      // riempio il vettore 'vector<TH2F*>& grafici_chi2'

            if (par.size() == 2) return;    //Se i parametri sono 2 è solo 1 il grafico da fare ed esco
        }
    }
}

//Grafici dei residui
void PlotGenerator::plot_residui(vector<TGraphErrors*>& grafici_residui, vector<double> par) {

    if (sigma_y.size() == 0) return;

    bool check = false;
    for (size_t g = 0; g < 3; g++)
    {
        if (check)
            return;         //se ho fatto il solo grafico dei residui in y poi non faccio altri grafici dei residui

        if (sigma_x.size() == 0) {
            g = 1;              //se non ci sono errori in x faccio il solo grafico dei residui in y
            check = true;
        }

        gStyle->SetOptStat(0);
        TGraphErrors* dati_plot = new TGraphErrors();

        string title = "";
        string asse_x = "";
        string asse_y = "";

        int dati_size = x.size();
        for (int i = 0; i < dati_size; i++) {

            double val_x = 0;
            double val_y = 0;
            double err_x = 0;
            double err_y = 0;

            double residuo_x = (sigma_x.size() != 0) ? x[i] - x_function(par, i) : 0;       //usa un metodo di bisezione per approssimare x(y)
            double residuo_y = y[i] - funzione_interpolante(par, x[i]);

            if (g == 0)              //residui in x
            {
                val_x = y[i];
                val_y = residuo_x;
                err_x = 0;
                err_y = (sigma_x.size() != 0) ? sigma_x[i] : 0;
                title = "Grafico Residui x";
                asse_x = "y";
                asse_y = "residui x";
            }
            else if (g == 1)         //residui in y
            {
                val_x = x[i];
                val_y = residuo_y;
                err_x = 0;
                err_y = (sigma_y.size() != 0) ? sigma_y[i] : 0;
                title = "Grafico Residui y";
                asse_x = "x";
                asse_y = "residui y";
                //cout << val_x << "\t" << val_y << "\t" << err_x << "\t" << err_y << endl;
            }
            else                     //residui in x e y
            {
                val_x = residuo_x;
                val_y = residuo_y;
                err_x = (sigma_x.size() != 0) ? sigma_x[i] : 0;
                err_y = (sigma_y.size() != 0) ? sigma_y[i] : 0;
                title = "Grafico Residui x-y";
                asse_x = "residui x";
                asse_y = "residui y";
            }

            dati_plot->SetPoint(i, val_x, val_y);
            dati_plot->SetPointError(i, err_x, err_y);
        }

        dati_plot->SetTitle(title.c_str());
        dati_plot->GetXaxis()->SetTitle(asse_x.c_str());
        dati_plot->GetYaxis()->SetTitle(asse_y.c_str());
        dati_plot->SetMarkerStyle(22);
        dati_plot->SetMarkerSize(0.7);
        dati_plot->SetLineWidth(1);

        grafici_residui.push_back(dati_plot);
    }

}


//Creazione dei canvas dei grafici delle distribuzioni del chi quadro
void PlotGenerator::canvas_chi_distribution(TCanvas*& c, TH2F*& grafico_chi2, vector<double> par, vector<double> sigma_par, int num_coppia) {

    c->cd();

    // Disegna come prima cosa in "sfondo" il grafico (poi sopra le linee)
    grafico_chi2->Draw("colz");

    //Defininizione dei margini
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.13);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);

    double i = num_coppia;
    double j = (num_coppia < 2) ? (num_coppia + 1) : 0;

    double x_i = par[i] - sigma_par[i];
    double x_f = par[i] + sigma_par[i];
    double y_i = par[j] - sigma_par[j];
    double y_f = par[j] + sigma_par[j];
    double x_min = grafico_chi2->GetXaxis()->GetXmin();
    double x_max = grafico_chi2->GetXaxis()->GetXmax();
    double y_min = grafico_chi2->GetYaxis()->GetXmin();
    double y_max = grafico_chi2->GetYaxis()->GetXmax();

    // Creazione delle linee di +- 1 sigma
    TLine* line1 = new TLine(x_i, y_min, x_i, y_max);
    TLine* line2 = new TLine(x_f, y_min, x_f, y_max);
    TLine* line3 = new TLine(x_min, y_i, x_max, y_i);
    TLine* line4 = new TLine(x_min, y_f, x_max, y_f);

    // Imposta lo stile delle linee
    line1->SetLineColor(kBlack);
    line2->SetLineColor(kBlack);
    line3->SetLineColor(kBlack);
    line4->SetLineColor(kBlack);

    line1->SetLineStyle(2);  // Linea tratteggiata
    line2->SetLineStyle(2);
    line3->SetLineStyle(2);
    line4->SetLineStyle(2);

    // Spessore delle linee
    line1->SetLineWidth(2); // Linea di spessore 2 pixel
    line2->SetLineWidth(2); // Linea di spessore 2 pixel
    line3->SetLineWidth(2); // Linea di spessore 2 pixel
    line4->SetLineWidth(2); // Linea di spessore 2 pixel

    // Disegna le linee sul Canvas
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");
    line4->Draw("same");



    c->Update();

}

//Creazione dei canvas dei grafici dei residui
void PlotGenerator::canvas_residui(TCanvas*& c, TGraphErrors*& grafico_residui, int num_grafico) {

    c->cd();

    // Disegna come prima cosa in "sfondo" il grafico (poi sopra le linee)
    grafico_residui->Draw("AP");

    //Defininizione dei margini
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);

    double x_min = grafico_residui->GetXaxis()->GetXmin();
    double x_max = grafico_residui->GetXaxis()->GetXmax();
    double y_min = grafico_residui->GetYaxis()->GetXmin();
    double y_max = grafico_residui->GetYaxis()->GetXmax();


    // Creazione delle linee di zero
    TLine* line_x = new TLine(x_min, 0.0, x_max, 0.0);
    TLine* line_y = new TLine(0.0, y_min, 0.0, y_max);

    // Imposta lo stile delle linee
    line_x->SetLineColor(kRed);
    line_y->SetLineColor(kRed);

    line_x->SetLineStyle(1);  // Linea continua
    line_y->SetLineStyle(1);

    // Spessore delle linee
    line_x->SetLineWidth(2); // Linea di spessore 2 pixel
    line_y->SetLineWidth(2); // Linea di spessore 2 pixel

    // Disegna le linee sul Canvas
    line_x->Draw("same");
    if (num_grafico == 2) line_y->Draw("same");

    c->Update();
}


// Funzione per validare un'estensione
bool PlotGenerator::isValidFormat(const std::string& extension) {
    // Lista delle estensioni supportate da ROOT
    static const std::set<std::string> validExtensions = {
        ".pdf", ".png", ".jpg", ".jpeg", ".eps", ".svg", ".root", ".C", ".gif", ".tiff"
    };

    // Controlla se l'estensione è nella lista di quelle valide
    return validExtensions.find(extension) != validExtensions.end();
}