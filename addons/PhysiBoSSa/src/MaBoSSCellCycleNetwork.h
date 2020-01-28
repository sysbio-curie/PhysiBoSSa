#ifndef __CELL_CYCLE_NETWORK_H
#define __CELL_CYCLE_NETWORK_H

#include "MaBoSSNetwork.h"

/**
 * \class CellCycleNetwork
 * \brief Cell Cycle directed by boolean network
 *
 *	Cell cycle and cell properties are controlled by the boolean network linked to the cell, calling MaBoSS instance.
 *	Takes as input current state of the cell (surrounding, acces to O2, contact with other cells, with ECM...) and read outputs from the netowrk (Migration, Proliferation...), and update cell state accordingly (cycle phase, migrative parameter...).
 *
 * Possible inputs and outputs depends on the structures of the network, try to put different possibilities.
 *
 *	\date 05/30/2017
 *	\author Gaelle Letort, Institut Curie
 */

class CellCycleNetwork
{
	private:
		/** \brief MaBoss Network doing the job */
		MaBoSSNetwork* maboss;
		
		/** \brief Vector of nodes state current value (0 or 1) */
		std::vector<bool> nodes;

		/** \brief time left before to update it again */
		double time_to_update;
	
		/** \brief choose a random update time, to asynchronize it between all cells 
		 *
		 * Set the next time at which to update the current cell's network. The time in between two udpates is chosen randomly in order to not update all cells together. */
		void set_time_to_update();

	public:
		/** \brief Constructor with given network instance */
		CellCycleNetwork(std::string bnd_file, std::string cfg_file);
		/** \brief Class destructor */
		~CellCycleNetwork();

		inline std::vector<bool>* get_nodes() {return &nodes;}
		inline double get_time_to_update() {return time_to_update;}
		inline int get_maboss_node_index( std::string name ) {return maboss->get_node_index(name);}
		inline void print_nodes() {maboss->print_nodes();}

		/** \brief Update MaboSS network states */
		void run_maboss();
};

#endif