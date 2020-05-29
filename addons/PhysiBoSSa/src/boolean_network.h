#ifndef __CELL_CYCLE_NETWORK_H
#define __CELL_CYCLE_NETWORK_H

#include "maboss_network.h"
#include "../../../core/PhysiCell_utilities.h"

/**
 *	\class BooleanNetwork
 *	\brief Boolean network handeler
 * 
 *	\details Boolean network class stores the values of a boolean network for a specific cell
 *
 *	\date 06/08/2020
 *	\author Gaelle Letort, Institut Curie
 *	\author Gerard Pradas, BSC-CNS
 */
class BooleanNetwork
{
	private:
		/** \brief MaBoss Network doing the job */
		MaBoSSNetwork maboss;
		
		/** \brief Vector of nodes state current value (0 or 1) */
		std::vector<bool> nodes;

		/** \brief Time left before to update it again */
		double time_to_update;
	
		/** \brief choose a random update time, to asynchronize it between all cells 
		 *
		 * Set the next time at which to update the current cell's network. The time in between two udpates is chosen randomly in order to not update all cells together. */
		inline void set_time_to_update(){this->time_to_update = (PhysiCell::UniformRandom()+0.5) * this->maboss.get_update_time_step();}

	public:
		/** 
		 * \brief Initialize a maboos network 
		 * \param bnd_file,cfg_file The configuration files from the boolean network
		 * \param time_step Time step between each MaBoSS simulation
		 */
		void initialize_boolean_network(std::string bnd_file, std::string cfg_file, double time_step);

		/** \brief Initialize a maboos network, with initial values and mutants */
		void initialize_boolean_network(std::string bnd_file, std::string cfg_file, double time_step, std::map<std::string, double> initial_values, std::map<std::string, double> mutations, std::map<std::string, double> parameters);
		
		/** \brief Reset nodes and time to update */
		void restart_nodes();

		/** \brief Update MaboSS network states */
		void run_maboss();
		
		/** 
		 * \brief Get nodes 
		 * \return The bool vector of the network values
		 */
		inline std::vector<bool>* get_nodes() {return &this->nodes;}
		
		/** 
		 * \brief Get time to update
		 * \return Time left to run next simulation 
		 */
		inline double get_time_to_update() {return this->time_to_update;}
		
		/**
		 * \brief Get value of a node by name
		 * \param name Node name existing in the boolean network
		 * \return Value of the node
		 */
		bool get_node_value( std::string name );

		/**
		 * \brief Set value of a node by name
		 * \param name Node name existing in the boolean network
		 * \param value Value to set the name node
		 * */
		void set_node_value( std::string name, bool value );

		/**
		 * \brief Get index of a node by name
		 * \param name Node name existing in the boolean network
		 * \return Index of the name node
		 */
		int get_node_index( std::string name );
		
		/** \brief Print name and values from the network */
		inline void print_nodes() {this->maboss.print_nodes(&this->nodes);}
		
		bool has_init() { return this->maboss.has_init(); }
};

#endif