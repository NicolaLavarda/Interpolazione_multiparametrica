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


extern vector<double> x, y, sigma_y;                             //Dati iniziali


//Funzione obbiettivo chi_quadro da usare nel calcolo delle derivate seconde per l'hessiana
double f(vector<double> parametri) {
    return f_chi_quadro(parametri);
}

// Funzione per il calcolo numerico della derivata prima tramite differenze finite con O(h^4)
double first_derivative_Oh4(function<double(vector<double>, vector<double>, int)> f, vector<double> params, int i) {
    double h = 1e-7;
    vector<double> p1 = params, p2 = params, p3 = params, p4 = params;

    p1[i] += 2 * h;  // f(p_i + 2h)
    p2[i] += h;      // f(p_i + h)
    p3[i] -= h;      // f(p_i - h)
    p4[i] -= 2 * h;  // f(p_i - 2h)

    return (-f(x, p1, 1) + 8 * f(x, p2, 1) - 8 * f(x, p3, 1) + f(x, p4, 1)) / (12 * h);
}


// Funzione per il calcolo numerico della derivata seconda tramite differenze finite
double second_derivative1(vector<double> params, int i, int j) {                       //function<double(const std::vector<double>&)> f
    double h = 1e-7;
    vector<double> p1 = params, p2 = params, p3 = params, p4 = params;

    p1[i] += h; p1[j] += h;  // f(p_i + h, p_j + h)
    p2[i] += h;              // f(p_i + h, p_j)
    p3[j] += h;              // f(p_i, p_j + h)
    // f(p_i, p_j) is simply f(params)
    return (f(p1) - f(p2) - f(p3) + f(params)) / (h * h);
}

// Calcolo derivate seconde miste
double second_derivative(vector<double> params, int i, int j) {
    double h = 1e-5;
    vector<double> p1 = params, p2 = params, p3 = params, p4 = params;

    p1[i] += h; p1[j] += h;  // f(p_i + h, p_j + h)
    p2[i] += h; p2[j] -= h;  // f(p_i + h, p_j - h)
    p3[i] -= h; p3[j] += h;  // f(p_i - h, p_j + h)
    p4[i] -= h; p4[j] -= h;  // f(p_i - h, p_j - h)

    return (f(p1) - f(p2) - f(p3) + f(p4)) / (4 * h * h);
}

// Calcolo derivate seconde pure
double second_derivative(vector<double> params, int i) {
    double h = 1e-5;
    vector<double> p1 = params, p2 = params, p3 = params;

    p1[i] += h;  // f(p_i + h)
    p2[i] -= h;  // f(p_i - h)

    return (f(p1) - 2 * f(params) + f(p2)) / (h * h);
}

// Funzione per calcolare l'intera matrice hessiana
vector<vector<double>> hessian(vector<double>& params) {
    int n = params.size();
    vector<vector<double>> H(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {  // Solo metà matrice (simmetrica)
            if (i != j) {
                H[i][j] = second_derivative(params, i, j);
                H[j][i] = H[i][j];
            }
            else
            {
                H[i][i] = second_derivative(params, i);
            }
        }
    }
    return H;
}






// Funzione per stampare una matrice
void stampaMatrice(const vector<vector<double>>& mat) {
    for (const auto& row : mat) {
        cout << "(   ";
        for (double value : row) {
            cout << left << setw(10) << setprecision(10) << value << "   ";
        }
        cout << ")" << endl;
    }
}

// Funzione per calcolare il determinante di una matrice
double determinante(const vector<vector<double>>& mat) {
    int n = mat.size();

    // Caso base: determinante di una matrice 1x1
    if (n == 1) return mat[0][0];

    double det = 0.0;
    for (int i = 0; i < n; i++) {
        // Creazione della sottomatrice per il cofattore
        vector<vector<double>> submat(n - 1, vector<double>(n - 1));

        for (int row = 1; row < n; row++) {
            int colIndex = 0;
            for (int col = 0; col < n; col++) {
                if (col == i) continue;
                submat[row - 1][colIndex++] = mat[row][col];
            }
        }

        // Ricorsione per calcolare il determinante della sottomatrice
        det += (i % 2 == 0 ? 1 : -1) * mat[0][i] * determinante(submat);
    }
    return det;
}

// Funzione per calcolare la matrice aggiunta (dei cofattori)
vector<vector<double>> matriceCofattori(const vector<vector<double>>& mat) {
    int n = mat.size();
    vector<vector<double>> cofmat(n, vector<double>(n));
    if (n == 1)
    {
        cofmat[0][0] = 1;
        return cofmat;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Creazione della sottomatrice per il cofattore
            vector<vector<double>> submat(n - 1, vector<double>(n - 1));

            int sub_i = 0;
            for (int row = 0; row < n; row++) {
                if (row == i) continue;
                int sub_j = 0;
                for (int col = 0; col < n; col++) {
                    if (col == j) continue;
                    submat[sub_i][sub_j++] = mat[row][col];
                }
                sub_i++;
            }

            // Calcolo del cofattore
            cofmat[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * determinante(submat);
        }
    }
    return cofmat;
}

// Funzione per trasporre una matrice
vector<vector<double>> trasposta(const vector<vector<double>>& mat) {
    int m = mat.size();             // Numero di righe della matrice originale
    int n = mat[0].size();           // Numero di colonne della matrice originale

    // Creazione della matrice trasposta di dimensioni n x m
    vector<vector<double>> tras(n, vector<double>(m));

    // Riempimento della matrice trasposta
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            tras[j][i] = mat[i][j];
        }
    }
    return tras;
}

// Funzione per calcolare l'inversa di una matrice
vector<vector<double>> inversa(const vector<vector<double>>& H) {
    int n = H.size();

    // Controllo che la matrice sia quadrata
    for (const auto& row : H) {
        if (row.size() != n) {
            throw invalid_argument("La matrice deve essere quadrata.");
        }
    }

    // Calcolo del determinante
    double det = determinante(H); //cout << det << endl;

    // Se il determinante è zero, la matrice non è invertibile
    if (det == 0) {
        throw invalid_argument("La matrice non è invertibile (determinante = 0).");
    }

    // Calcolo della matrice dei cofattori
    vector<vector<double>> cofmat = matriceCofattori(H);
    //cout << "Matrice cofattoti:" << endl;
    //stampaMatrice(cofmat);

    // Trasporre la matrice dei cofattori (matrice aggiunta)
    vector<vector<double>> adj = trasposta(cofmat);
    //cout << "Matrice cofattoti trasposta:" << endl;
    //stampaMatrice(adj);

    // Calcolo della matrice inversa (aggiunta / determinante)
    vector<vector<double>> invmat(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            invmat[i][j] = adj[i][j] / det;
        }
    }

    return invmat;
}
