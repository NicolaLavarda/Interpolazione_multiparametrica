#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <functional>

using namespace std;

double f(vector<double> parametri);

double first_derivative_Oh4(function<double(vector<double>, vector<double>, int)> f, vector<double> params, int i);

double second_derivative1(vector<double> params, int i, int j);

double second_derivative(vector<double> params, int i, int j);

double second_derivative(vector<double> params, int i);

vector<vector<double>> hessian(vector<double>& params);




void stampaMatrice(const vector<vector<double>>& mat);

double determinante(const vector<vector<double>>& mat);

vector<vector<double>> matriceCofattori(const vector<vector<double>>& mat);

vector<vector<double>> trasposta(const vector<vector<double>>& mat);

vector<vector<double>> inversa(const vector<vector<double>>& H);


#endif