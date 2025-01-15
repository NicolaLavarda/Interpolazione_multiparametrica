#define ChiSquareMinimizer_h
#ifdef ChiSquareMinimizer_h

#include "covering.h"

#include <vector>
#include <string>
#include <map>

class ChiSquareMinimizer {
public:
	ChiSquareMinimizer(std::map<std::string, bool>& options, int num_val);

	void begin(std::vector<double>& par_best);

	void compute(std::vector<double>& par_best);

	void end(std::vector<double>& par_best);



	void thread1(std::vector<double>& par_best);
	void thread2(std::vector<double>& par_best);


private:

	bool faster, approx, complex, retta;	  // opzioni da riga di comando

	std::vector<double> par_best;		 //parametri migliori in assoluto
	double chi_quadro_min = 1e30;        //Chi quadro minimo assoluto (in ogni istante di tutto il programma)

	bool errore_lin = false;		// serve per 'compute'
	int cicle_programms = 0;        // Numero di cicli di riesecuzione del programma

	int x_size;

	//bool finish = false;

	covering* c_generator = nullptr;

};

#endif
