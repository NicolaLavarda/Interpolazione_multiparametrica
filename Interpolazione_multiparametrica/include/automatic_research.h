#ifndef AUTOMATIC_RESEARCH_H
#define AUTOMATIC_RESEARCH_H

#include "interpolator.h"

#include <vector>

class AutomaticResearch {
public:
	AutomaticResearch(std::vector<double> par, bool output);

	void beginJob();

	void endJob(std::vector<double>& par);

private:

	void compute(std::vector<double>& par_auto, int n);


	std::vector<double> par;	// parametri da trovare/migliorare
	bool output;	// flag per stampare a schermo calcoli intermedi

	// Riferimento all'istanza Singleton
	Interpolator& i_generator = Interpolator::getInstance();

	// Vettore che contiene i vettori di parametri di volta in volta migliorati
	std::vector<std::vector<double>> par_improved;

	// indice del primo parametro automatico
	int first;

};


#endif