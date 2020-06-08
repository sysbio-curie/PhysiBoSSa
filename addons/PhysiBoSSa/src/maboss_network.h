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
#include "../../../core/PhysiCell_utilities.h"

class MaBoSSNetwork
{
	private:
		/** \brief MaBoSS instances: network */
		Network* network;
		/** \brief MaBoSS instances: configurations */
		RunConfig* config;

		StochasticSimulationEngine* engine;
		NetworkState state;
		NetworkState output_mask;
		
		/** \brief Time step to update the cycle */
		double update_time_step = 12.0;
		
		/** \brief Real time to update, after applying noise */
		double time_to_update;

		double scaling = 1.0;
		
		/** \brief Initial value probabilities, by node */
		std::map< std::string, double > initial_values;
		
		/** \brief Mutations to apply to the network */
		std::map< std::string, double > mutations;
	
		std::map< std::string, Node*> nodesByName;
		std::map< std::string, const Symbol*> parametersByName;
	
		inline void set_time_to_update(){this->time_to_update = (PhysiCell::UniformRandom()+0.5) * this->get_update_time_step();}

	
	public:
	
		/** Constructor */
		MaBoSSNetwork() {
			network = NULL;
			config = NULL;
			engine = NULL;
			this->nodesByName.clear();
		}
		
		/** Desctructor */
		~MaBoSSNetwork() {
			delete this->engine;
			this->engine = NULL;
			delete this->network;
			this->network = NULL;
			delete this->config;
			this->config = NULL;
		}
		
		/** \brief Initialize network */
		void init_maboss( std::string networkFile, std::string configFile);

		bool has_init() const { return network != NULL && config != NULL; }
		void mutate(std::map<std::string, double> mutations);

		void set_initial_values(std::map<std::string, double> initial_values)
		{ this->initial_values = initial_values; }

		void set_parameters(std::map<std::string, double> parameters);

		double get_parameter_value(std::string name);
		void set_parameter_value(std::string name, double value);
		
		/** \brief Restart a vector of bools, to the init values of the network */
		void restart_node_values();

		/** \brief Run the current network*/
		void run_simulation();

		/** \brief Return node of given name current value
		 *
		 * Return -1 if node doesn't exit \n
		 * Return 0 if node is 0 \n
		 * Return 1 if node is 1 */
		// int get_maboss_node_index( std::string name );

		bool has_node( std::string name );
		void set_node_value(std::string name, bool value);
		bool get_node_value(std::string name);
		std::string get_state();
		
		/** \brief Return update time value */
		inline double get_update_time_step(){ return this->update_time_step; }

		/** \brief Setter update time */
		inline void set_update_time_step(double time_step) { this->update_time_step = time_step;}
		
		
		/** \brief Get time to update*/
		inline double get_time_to_update() {return this->time_to_update;}
		
		/** \brief Change simulation mode */
		inline void set_discrete_time(bool discrete_time, double time_tick) { 
			this->engine->setDiscreteTime(discrete_time); this->engine->setTimeTick(time_tick); 
		}

		inline void set_scaling(double scaling) { this->scaling = scaling; }
		
		/** \brief Print current state of all the nodes of the network */
		void print_nodes();
		
		

};

#endif
