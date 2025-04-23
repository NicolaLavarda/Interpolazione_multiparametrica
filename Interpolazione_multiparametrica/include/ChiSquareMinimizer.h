#define ChiSquareMinimizer_h
#ifdef ChiSquareMinimizer_h

#include "covering.h"

#include <vector>
#include <string>
#include <map>

class ChiSquareMinimizer {
public:
	ChiSquareMinimizer(std::map<std::string, bool>& options);

	void begin(std::vector<double>& par_best);

	void compute(std::vector<double>& par_best);

	void end(std::vector<double>& par_best);

private:

	void decreseSteps(const std::vector<double>& par_best, const std::vector<double>& par_best_prec);

	void computeCovering(std::vector<double>& par_best, int cicle_programm);

	bool faster, approx, complex;	  // opzioni da riga di comando

	//std::vector<double> par_best;		 //parametri migliori in assoluto
	double chi_quadro_min = 1e30;        //Chi quadro minimo assoluto (in ogni istante di tutto il programma)

	// variabili per 'computeCovering'
	covering* c_generator = nullptr;
	std::vector<double> steps;
	bool update_steps = false;

};

#endif
