#ifndef PLOT_H
#define PLOT_H

#include <vector>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TCanvas.h>


std::vector<std::vector<double>> calcolo_matrice_chi_quadri(std::vector<double> par, std::vector<double> sigma_par, double num_sigma, double punti_per_parametro, int primo_par, int secondo_par);


void plot_function(TMultiGraph*& grafico_dati_interpolazione, std::vector<double> par, std::vector<double> sigma_par);


void plot_chi_distribution(std::vector<TH2F*>& grafici_chi2, std::vector<double> par, std::vector<double> sigma_par);

void plot_residui(std::vector<TGraphErrors*>& grafici_residui, std::vector<double> par);

void canvas_chi_distribution(TCanvas*& c, TH2F*& grafico_chi2, std::vector<double> par, std::vector<double> sigma_par, int num_coppia);

void canvas_residui(TCanvas*& c, TGraphErrors*& grafico_residui, int num_grafico);

#endif