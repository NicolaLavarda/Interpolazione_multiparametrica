#ifndef INPUT_H
#define INPUT_H

#include "interpolator.h"

#include <vector>
#include <string>
#include <map>

class input {
public:
	input(int argc, char* argv[]);

	void compute(std::vector<double>& par, int& num_a, std::map<std::string, bool>& options,
				  std::vector<double>& x, std::vector<double>& sigma_x,
				  std::vector<double>& y, std::vector<double>& sigma_y);

	bool improved(const std::string filePath, std::vector<double> par_best);

	std::vector<double> GetParametersFromFile();

	std::vector<double> GetErrorsFromFile();



private:
	// restituisce il valore di un numero contenuto nel file 'filePath', nella riga in cui si trova 'Target', dopo il carattere 'AfterThat'
	double GetValueFromFile(const std::string filePath, const std::string target, const std::string afterThat);

	double GetParameter(const std::string parameter);

	double GetErrorOfParameter(const std::string parameter);


	// Prendo l'istanza al Singleton 'Interpolator'
	Interpolator& i_generator = Interpolator::getInstance();

	std::string filePath;

	// Argomenti di 'input'
	int argc;
	char** argv;

	// Argomenti di 'beginJob'
	std::vector<double> x, sigma_x, y, sigma_y;
	std::vector<double> par;


	std::vector<std::string> name_parameters = { "a","b","c", "d", "e" };
};





#endif