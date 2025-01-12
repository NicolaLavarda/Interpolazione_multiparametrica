#ifndef PLOT_H
#define PLOT_H

#include "interpolator.h"

#include <vector>
#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>


class PlotGenerator {
public:

    PlotGenerator(std::vector<double> par, std::vector<double> x, std::vector<double> sigma_x, std::vector<double> y, std::vector<double> sigma_y, std::string file_name);

    void compute_plot_function();

    void compute_plot_chi_distribution();

    void compute_plot_Residuals();

    void save(std::string format_file);

private:

    // Generation of Graphs

    void plot_function(TMultiGraph*& grafico_dati_interpolazione, std::vector<double> par);

    void plot_chi_distribution(std::vector<TH2F*>& grafici_chi2, std::vector<double> par, std::vector<double> sigma_par);

    void plot_residui(std::vector<TGraphErrors*>& grafici_residui, std::vector<double> par);


    // Update Canvas

    void canvas_plot_function(TCanvas*& c, TLegend*& legend, TMultiGraph*& grafico);

    void canvas_chi_distribution(TCanvas*& c, TLegend*& legend, TH2F*& grafico_chi2, std::vector<double> par, std::vector<double> sigma_par, int num_coppia);

    void canvas_residui(TCanvas*& c, TLegend*& legend, TGraphErrors*& grafico_residui, int num_grafico);


    // Utilities

    std::vector<std::vector<double>> calcolo_matrice_chi_quadri(std::vector<double> par, std::vector<double> sigma_par, double num_sigma, double punti_per_parametro, int primo_par, int secondo_par);

    bool isValidFormat(const std::string& extension);


    // Dati inizializzati dal costruttore
    std::vector<double> par;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> sigma_x;
    std::vector<double> sigma_y;

    // Grafici da plottare
    TMultiGraph* grafico_dati_interpolazione;
    std::vector<TH2F*> grafici_chi2;
    std::vector<TGraphErrors*> grafici_residui;

    // Vettore di tutte le TLegend da associare ai relativi Canvas
    std::vector<TLegend*> Legend_grafici;

    // Vettore che contiene i puntatori a tutti i Canvas che mi servono
    std::vector<TCanvas*> Canvas_grafici;

    //formato del file in cui salvare i grafici
    std::string format = ".png";

    //Nome base di tutti i grafici
    std::string base_name;

    // Riferimento all'istanza Singleton
    Interpolator& i_generator = Interpolator::getInstance();

};



#endif