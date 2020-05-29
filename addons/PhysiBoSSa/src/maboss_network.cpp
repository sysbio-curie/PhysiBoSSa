#include "maboss_network.h"
#include "../../../core/PhysiCell_utilities.h"

/* Default constructor */
void MaBoSSNetwork::init_maboss( std::string networkFile, std::string configFile)
{
	if (this->network != NULL) {
		delete this->network;
	}
	
	if (this->config != NULL) {
		delete this->config;
	}
	
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->parse(this->network, configFile.c_str());

	IStateGroup::checkAndComplete(this->network);

	engine = new StochasticSimulationEngine(this->network, this->config, PhysiCell::UniformInt());

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
	delete this->engine;
	this->engine = NULL;
	delete this->network;
	this->network = NULL;
	delete this->config;
	this->config = NULL;
}

void MaBoSSNetwork::mutate(std::map<std::string, double> mutations) 
{
	for (auto mutation : mutations) {
		Node * node = this->network->getNode(mutation.first);
		node->mutate(mutation.second);
	}
}

void MaBoSSNetwork::set_parameters(std::map<std::string, double> parameters) 
{	
	SymbolTable* symbol_table = this->network->getSymbolTable();
	
	for (auto parameter: parameters) {
		const Symbol* symbol_paraneter = symbol_table->getSymbol(parameter.first);
		symbol_table->setSymbolValue(symbol_paraneter, parameter.second);
	}
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
	this->network->initStates(network_state, engine->random_generator);

	int i = 0;
	std::vector<Node *> nodes = this->network->getNodes();
	(*node_values).resize(nodes.size());
	for (auto node : nodes)
	{
		(*node_values)[i] =  network_state.getNodeState( node ) ;;
		i++;
	}
	
	for (auto initial_value : initial_values) {
		(*node_values)[node_names[initial_value.first]] = PhysiCell::UniformRandom() < initial_value.second;
	}
}

/* Run the current network */
void MaBoSSNetwork::run_simulation(std::vector<bool>* node_values)
{	
	NetworkState_Impl state = this->load_state(node_values);

	state = engine->run(&state, NULL);

	this->recover_state(state, node_values);
}

/* Return the index of node based on node's name */
int MaBoSSNetwork::get_maboss_node_index( std::string name )
{
	auto res = this->node_names.find(name);
	if ( res != this->node_names.end() )
		return res->second;
		
	std::string err_msg = "A node with name " + name + " does not exist in the network.";
	throw std::invalid_argument(err_msg);
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