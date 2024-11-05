#ifndef BRANCHING_COVER
#define BRANCHING_COVER

#include "interpolator.h"

#include <vector>
#include <string>
#include <map>

class BranchingCover {
public:
	BranchingCover(std::vector<double>& par_best);

	void doBetter(std::vector<double>& par_best);

	~BranchingCover();

private:

	void doBetter(std::vector<double>& par_best, int address);

	int register_par_address(int i, int address);

	bool up_down(std::vector<double>& par_i, std::vector<double> grad);

	std::vector<std::vector<double>> division_direction(const std::vector<double> par);

	std::vector<double> grad_f_chi_quadro(const std::vector<double> par);


	// Accesso ad una nuova istanza già settata
	Interpolator* i_generator = Interpolator::getNewInstance();

	int par_size;

	// contenitore dei parametri migliorati nelle diverse direzioni/divisioni
	std::vector<std::vector<double>> par_div;

	// mappa dove associare alla posizione di un vettore nella ramificazione l'indice corrispondente in 'par_div'
	std::map<int, int> address_par;

	// livello corrente profondità ramificazione (numero di chiamate di 'doBetter')
	int level;

	// numero livelli massimi
	int max_level;


};


#endif