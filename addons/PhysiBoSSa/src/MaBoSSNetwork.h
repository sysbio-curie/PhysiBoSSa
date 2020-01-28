#ifndef _MaBoSS_Net_h_
#define _MaBoSS_Net_h_

/**
 * \brief Interface with MaBoSS software
 *
 * One network object contains network configs.
 * Handle initialization, running it...
 *
 * Created on 05/19/2017
 * G. Letort, Institute Curie
 */

#include "StochasticSimulationEngine.h"
#include "BooleanNetwork.h"
#include "RunConfig.h"
#include "utils.h"

class MaBoSSNetwork
{
	private:
		/** \brief MaBoSS instances: network */
		Network* network;
		/** \brief MaBoSS instances: configurations */
		RunConfig* config;
		/** \brief MaBoSS instances: state */
		NetworkState_Impl* state;

		/** \brief Time step to update the cycle */
		double update_time_step;

		/** \brief Names and indices of network nodes */
		std::map< std::string, int > node_names;

	public:
		/** \brief Constructor */
		MaBoSSNetwork( std::string networkFile, std::string configFile );
		/** \brief Destructor */
		~MaBoSSNetwork();
		
		inline int nodes_lenght() { return node_names.size(); }
		/** \brief Return update time value */
		inline double get_update_time_step(){ return update_time_step; }

		void load_state(std::vector<bool>* input);
		void recover_state(std::vector<bool>* output);

		/** \brief Run the current network*/
		void run();

		void run(std::vector<bool>* node_values);
		
		/** \brief Print current state of all the nodes of the network */
		void print_nodes();

		/** \brief Return node of given name current value
		 *
		 * Return -1 if node doesn't exit \n
		 * Return 0 if node is 0 \n
		 * Return 1 if node is 1 */
		int get_node_index( std::string name );
};

#endif