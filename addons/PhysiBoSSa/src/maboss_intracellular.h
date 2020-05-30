#ifndef _MaBoSS_Intracellular_h_
#define _MaBoSS_Intracellular_h_

#include <string>
#include <map>
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "boolean_network.h"

class MaBoSSIntracellular : public PhysiCell::Intracellular {
 private:
 public:
	
	std::string bnd_filename;
	std::string cfg_filename;
	
	double time_step = 1000.;
	bool discrete_time = false;
	double time_tick = 0.5;
	
	std::map<std::string, double> initial_values;
	std::map<std::string, double> mutations;
	std::map<std::string, double> parameters;
	
	BooleanNetwork network;
	
	double next_physiboss_run = 0;

	MaBoSSIntracellular();
	
	MaBoSSIntracellular(pugi::xml_node& node);
	
	MaBoSSIntracellular(MaBoSSIntracellular* copy);
	
	
	Intracellular* getIntracellularModel() {
		return static_cast<Intracellular*>(this);
	}
	
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	void start() {
		this->network.restart_nodes();
	}
	
	void update() {
		this->network.run_maboss();
		this->next_physiboss_run += this->network.get_time_to_update();
	}
	
	bool need_update() {
		return PhysiCell::PhysiCell_globals.current_time >= this->next_physiboss_run;
	}
};

MaBoSSIntracellular* getMaBoSSModel(PhysiCell::Phenotype& phenotype);

#endif