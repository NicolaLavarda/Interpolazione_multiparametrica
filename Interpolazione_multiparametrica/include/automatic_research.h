#ifndef AUTOMATIC_RESEARCH_H
#define AUTOMATIC_RESEARCH_H

#include "interpolator.h"

#include <vector>
#include <string>
#include <memory>

class AutomaticResearch {
public:
	AutomaticResearch(std::vector<double> par, bool output);

	void beginJob();

	void doBetter();

	void endJob(std::vector<double>& par);

	~AutomaticResearch();

private:

	void compute(std::vector<double>& par_auto, int n);

	void print(const std::string name, std::vector<double> par);


	std::vector<double> par;	// parametri da trovare/migliorare
	bool output;	// flag per stampare a schermo calcoli intermedi

	// Accesso ad una nuova istanza già settata
	Interpolator* i_generator = Interpolator::getNewInstance();

	// Vettore che contiene i vettori di parametri di volta in volta migliorati
	std::vector<std::vector<double>> par_improved;

	// indice del primo parametro automatico
	int first;

	//range e valori della ricerca automatica
	double chi_min_auto = 1e25;
	double min = -1e10;
	double max = 1e10;
	double min_ordine = 0.001;       //dev'essere positivo (>0)
	int par_size;

};


#endif