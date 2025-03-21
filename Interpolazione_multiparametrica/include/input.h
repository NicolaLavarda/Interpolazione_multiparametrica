#ifndef INPUT_H
#define INPUT_H

#include "interpolator.h"

#include <vector>
#include <string>
#include <map>

class input {
public:

	// Dati iniziali
	struct Data {
		std::vector<double> x;
		std::vector<double> sigma_x;
		std::vector<double> y;
		std::vector<double> sigma_y;
	};

	struct Interpolation {
		std::vector<double> par;
		std::string interpolating_function;
		std::string filePath;
		std::vector<std::string> name_parameters = { "a", "b", "c", "d", "e" };
	};


	// Costruttore per utilizzare info da riga di comando
	input(int argc, char* argv[]);

	// Metodi per riempire le due structs e 'options' con info da riga di comando ('argc' e 'argv')
	void compute(std::map<std::string, bool>& options);

	// Metodi per salvare dati in structs
	void SetData(std::vector<double>& x, std::vector<double>& sigma_x,
				 std::vector<double>& y, std::vector<double>& sigma_y);
	void SetInterpolation(std::vector<double> par, std::string interpolating_function, std::string filePath);

	// Metodi per ottenere le structs con i dati
	Data GetData();
	Interpolation GetInterpolation();

	bool improved(const std::string filePath, std::vector<double> par_best);

	std::vector<double> GetParametersFromFile();

	std::vector<double> GetErrorsFromFile();

	void help();

private:
	// restituisce il valore di un numero contenuto nel file 'filePath', nella riga in cui si trova 'Target', dopo il carattere 'AfterThat'
	double GetValueFromFile(const std::string filePath, const std::string target, const std::string afterThat);

	double GetParameter(const std::string parameter);

	double GetErrorOfParameter(const std::string parameter);

	// Argomenti di 'input'
	int argc;
	char** argv;

	// struct
	Data data;
	Interpolation interpolation;

};


#endif