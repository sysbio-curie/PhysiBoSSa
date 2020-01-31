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

		/** \brief Time step to update the cycle */
		double update_time_step = 12;

		/** \brief Names and indices of network nodes */
		std::map< std::string, int > node_names;

		/** \brief Given a vectors of bools, returns the corresponding state */
		NetworkState_Impl load_state(std::vector<bool>* input);

		/** \brief Given a state and a vectors of bools, writes the states into the vector */
		void recover_state(NetworkState_Impl state, std::vector<bool>* output);

	public:
		/** \brief Initialize network */
		void init_maboss( std::string networkFile, std::string configFile);

		/** \brief Destructor */
		void delete_maboss();

		/** \brief Restart a vector of bools, to the init values of the network */
		void restart_node_values(std::vector<bool>* node_values);

		/** \brief Run the current network*/
		void run_simulation(std::vector<bool>* node_values);

		/** \brief Return node of given name current value
		 *
		 * Return -1 if node doesn't exit \n
		 * Return 0 if node is 0 \n
		 * Return 1 if node is 1 */
		int get_node_index( std::string name );

		/** \brief Return update time value */
		inline double get_update_time_step(){ return update_time_step; }
		
		/** \brief Print current state of all the nodes of the network */
		void print_nodes(std::vector<bool>* node_values);
};

#endif
