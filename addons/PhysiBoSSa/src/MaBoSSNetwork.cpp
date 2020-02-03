#include "MaBoSSNetwork.h"

/* Default constructor */
void MaBoSSNetwork::init_maboss( std::string networkFile, std::string configFile)
{
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->parse(this->network, configFile.c_str());

	IStateGroup::checkAndComplete(this->network);

	int i = 0;
	std::vector<Node *> nodes = this->network->getNodes();
	for (auto node : nodes)
	{
		this->node_names[ node->getLabel() ] = i;
		i++;
	}
}

/** Default estructor */
void MaBoSSNetwork::delete_maboss()
{
	delete this->network;
	this->network = NULL;
	delete this->config;
	this->config = NULL;
}

/* Creates a NetworkState_Impl from the input */
NetworkState_Impl MaBoSSNetwork::load_state(std::vector<bool>* input)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	NetworkState state_to_load;
	for (auto node : nodes)
	{
		state_to_load.setNodeState(node, (NodeState) (*input)[i]);
		i ++;
	}
	return state_to_load.getState();
}

/* Transfer state values to a vector of bools */
void MaBoSSNetwork::recover_state(NetworkState_Impl state, std::vector<bool>* output)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	NetworkState state_to_retrieve = (NetworkState) state;
	for ( auto node: nodes )
	{
		(*output)[i] = state_to_retrieve.getNodeState( node ) ;
		i++;
	}
}

/* Reset a vector of bools to the init state of the network */
void MaBoSSNetwork::restart_node_values(std::vector<bool>* node_values)
{
	NetworkState network_state;
	RandomGeneratorFactory *randgen_factory = this->config->getRandomGeneratorFactory();
  	RandomGenerator *random_generator = randgen_factory->generateRandomGenerator(UniformInt());
	this->network->initStates(network_state, random_generator);

	int i = 0;
	std::vector<Node *> nodes = this->network->getNodes();
	(*node_values).resize(nodes.size());
	for (auto node : nodes)
	{
		(*node_values)[i] =  network_state.getNodeState( node ) ;;
		i++;
	}
}

/* Run the current network */
void MaBoSSNetwork::run_simulation(std::vector<bool>* node_values)
{	
	NetworkState_Impl state = this->load_state(node_values);

	StochasticSimulationEngine* engine = new StochasticSimulationEngine(this->network, this->config);
	engine->setSeed(UniformInt());
	state = engine->run(&state, NULL);
	delete engine;

	this->recover_state(state, node_values);
}

/* Return the index of node based on node's name */
int MaBoSSNetwork::get_node_index( std::string name )
{
	auto res = this->node_names.find(name);
	if ( res != this->node_names.end() )
		return res->second;
	return -1;	
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes(std::vector<bool>* node_values)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << (*node_values)[i] << "; ";
		i++;
	}
	std::cout << std::endl;
}