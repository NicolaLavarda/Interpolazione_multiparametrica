#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <functional>
#include <fstream>

using namespace std;

double f(const vector<double>& parametri);

double second_derivative(const vector<double>& params, int i, int j);

double second_derivative(const vector<double>& params, int i);

vector<vector<double>> hessian(const vector<double>& params);




void stampaMatrice(const vector<vector<double>>& mat, std::ostream& output);

double determinante(const vector<vector<double>>& mat);

vector<vector<double>> matriceCofattori(const vector<vector<double>>& mat);

vector<vector<double>> trasposta(const vector<vector<double>>& mat);

vector<vector<double>> inversa(const vector<vector<double>>& H);


#endif