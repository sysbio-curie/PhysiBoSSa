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
	double time_step;
	std::map<std::string, double> initial_values;
	std::map<std::string, double> mutations;
	std::map<std::string, double> parameters;
	
	BooleanNetwork network;

	MaBoSSIntracellular();
	
	MaBoSSIntracellular(pugi::xml_node& node);
	
	MaBoSSIntracellular(MaBoSSIntracellular* copy);
	
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	Intracellular* getIntracellularModel() {
		return static_cast<Intracellular*>(this);
	}
};

MaBoSSIntracellular* getMaBoSSModel(PhysiCell::Phenotype& phenotype);

#endif