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

	void computeCovering(std::vector<double>& par_best, int cicle_programms);

	bool faster, approx, complex, retta;	  // opzioni da riga di comando

	std::vector<double> par_best;		 //parametri migliori in assoluto
	double chi_quadro_min = 1e30;        //Chi quadro minimo assoluto (in ogni istante di tutto il programma)

	bool errore_lin = false;		// serve per 'compute'
	int cicle_programms = 0;        // Numero di cicli di riesecuzione del programma

	covering* c_generator = nullptr;

};

#endif
