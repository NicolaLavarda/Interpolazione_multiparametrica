#ifndef BRANCHING_COVER
#define BRANCHING_COVER

#include "interpolator.h"

#include <vector>
#include <string>
#include <map>

class BranchingCover {
public:
	BranchingCover(std::vector<double>& par_best);

	void compute(std::vector<double>& par_best);

	~BranchingCover();

private:

	void doBetter(std::vector<double>& par_best, int address);

	int register_par_address(int i, int address);

	bool go_on(std::vector<double>& par_i, std::vector<double> div);

	bool branch_kill(int address);

	int mother_address(int son_address);

	int GetLevel(int address);

	std::vector<std::vector<double>> division_direction(int address);


	// Accesso ad una nuova istanza già settata
	Interpolator* i_generator = Interpolator::getNewInstance();

	int par_size;

	// contenitore dei parametri migliorati nelle diverse direzioni/divisioni		(es. (1 2 3.5) )
	std::vector<std::vector<double>> par_div;

	// contenitore dei passi fatti per trovare i rispettivi parametri 'par_div'		(es. (0 0 0.5) )
	std::vector<std::vector<double>> div_steps;

	// contenitore (ordinato secondo 'address_par') dei parametri che sono migliorati (True) o meno (False) rispetto ai parametri madre che li hanno generati
	std::vector<bool> improvements;

	// mappa dove associare alla posizione di un vettore nella ramificazione l'indice corrispondente in 'par_div'
	std::map<int, int> address_par;

	// livello corrente profondità ramificazione (numero di chiamate di 'doBetter')
	int level;

	// numero livelli massimi
	int max_level;

	// numero di ramificazioni che falliscono dopo le quali muore il ramo (non si ridivide)
	double death_rate;

	// passi iniziali
	std::vector<std::vector<double>> div;

	// passo ramificazioni iniziali (es. 0.2 = 20% del parametro)
	double step;


};


#endif