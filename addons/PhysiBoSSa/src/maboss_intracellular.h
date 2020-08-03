#ifndef _MaBoSS_Intracellular_h_
#define _MaBoSS_Intracellular_h_

#include <string>
#include <map>
#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "maboss_network.h"

/**
 *	\class MaBoSSIntracellular
 *	\brief Interface with PhysiCell
 * 
 *	\details The MaBoSSIntracellular is an implementation of the PhysiCell::Intracellular abstract class
 *           In handles the XML parsing, object copying, results export, update checking and running simulations
 *
 *	\date 03/08/2020
 *	\author Gaelle Letort, Institut Curie
 *	\author Gerard Pradas, BSC-CNS
 *	\author Vincent Noel, Institut Curie
 */

class MaBoSSIntracellular : public PhysiCell::Intracellular {
 private:
 public:
	
	static long counter;
	
	std::string bnd_filename;
	std::string cfg_filename;
	
	double time_step = 12;
	bool discrete_time = false;
	double time_tick = 0.5;
	double scaling = 1.0;
	
	std::map<std::string, double> initial_values;
	std::map<std::string, double> mutations;
	std::map<std::string, double> parameters;
	
	MaBoSSNetwork maboss;

	double next_physiboss_run = 0;

	MaBoSSIntracellular();
	
	MaBoSSIntracellular(pugi::xml_node& node);
	
	MaBoSSIntracellular(MaBoSSIntracellular* copy);
	
	Intracellular* clone() {
		return static_cast<Intracellular*>(new MaBoSSIntracellular(this));
	}
	Intracellular* getIntracellularModel() {
		return static_cast<Intracellular*>(this);
	}
	
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	/** 
	 * \brief Initialize the intracellular model
	 */
	void start() {
		this->maboss.restart_node_values();
	}
	
	/** 
	 * \brief Update the intracellular model
	 */
	void update() {
		this->maboss.run_simulation();
		this->next_physiboss_run += this->maboss.get_time_to_update();
	}
	
	/** 
	 * \brief Check if the intracellular model needs to be updated
	 */
	bool need_update() {
		return PhysiCell::PhysiCell_globals.current_time >= this->next_physiboss_run;
	}
	
	/** 
	 * \brief Check is a node exists
	 */
	bool has_node(std::string name) {
		return this->maboss.has_node(name);
	}
	
	/** 
	 * \brief Returns the boolean value of a node
	 */
	bool get_boolean_node_value(std::string name) {
		return this->maboss.get_node_value(name);
	}
	
	/** 
	 * \brief Sets the boolean value of a node
	 */
	void set_boolean_node_value(std::string name, bool value) {
		this->maboss.set_node_value(name, value);
	}
	
	/** 
	 * \brief Return the value of a model parameter
	 */
	double get_parameter_value(std::string name) {
		return this->maboss.get_parameter_value(name);
	}
	
	/** 
	 * \brief Set the value of a model parameter
	 */
	void set_parameter_value(std::string name, double value) {
		this->maboss.set_parameter_value(name, value);
	}
	
	/** 
	 * \brief Returns the state of the intracellular model as a string
	 */
	std::string get_state() {
		return this->maboss.get_state();
	}
	
	/** 
	 * \brief Save physiboss state of this intracellular in csv
	 */
	static void save_PhysiBoSS(std::string path, std::string index);

	/** 
	 * \brief Prints the current value of the intracellular model
	 */
	void print_current_nodes(){
		this->maboss.print_nodes();
	};
};

#endif