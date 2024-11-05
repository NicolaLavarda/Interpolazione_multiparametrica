#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>

using namespace std;

//defined in main.cpp
double funzione_interpolante(vector<double> x, vector<double> par, int i);

#endif













/*
#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>
#include <string>
#include "exprtk.hpp"

using namespace std;

class FunzioneCompilata {
public:
    FunzioneCompilata(const std::string& expression_str);
    double operator()(const std::vector<double>& x, const std::vector<double>& par, int i);

private:
    double x_;
    std::vector<double> par_{ 3 };

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;
};

void initFunzioneCompilata(const std::string& expression_str);
double funzione_interpolante(const std::vector<double>& x, const std::vector<double>& par, int i);

#endif
*/









/*
#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <vector>
#include <string>
#include "exprtk.hpp"

using namespace std;

class FunzioneCompilata {
public:
    FunzioneCompilata(const std::string& expression_str);
    double operator()(const std::vector<double>& x, const std::vector<double>& par, int i);

private:
    double x_;
    std::vector<double> par_{ 3 };

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;
};

// Definizione globale di FunzioneCompilata
extern FunzioneCompilata funzione;

#endif
*/