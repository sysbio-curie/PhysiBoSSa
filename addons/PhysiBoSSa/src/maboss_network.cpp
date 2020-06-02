#include "maboss_network.h"

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

	this->update_time_step = this->config->getMaxTime();
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

/* Reset a vector of bools to the init state of the network */
void MaBoSSNetwork::restart_node_values()
{
	// NetworkState network_state;
	this->network->initStates(state, engine->random_generator);
	
	for (auto initial_value : initial_values) {
		state.setNodeState(network->getNode(initial_value.first), PhysiCell::UniformRandom() < initial_value.second);
	}
	
	this->set_time_to_update();
}

/* Run the current network */
void MaBoSSNetwork::run_simulation()
{	
	engine->setMaxTime(time_to_update/scaling);
	state = engine->run(state, NULL);
	this->set_time_to_update();

}

bool MaBoSSNetwork::has_node( std::string name ) {
	return network->isNodeDefined(name);
}

void MaBoSSNetwork::set_node_value(std::string name, bool value) {
	state.setNodeState(network->getNode(name), value);
}

bool MaBoSSNetwork::get_node_value(std::string name) {
	return state.getNodeState(network->getNode(name));
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes()
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; ";
		i++;
	}
	std::cout << std::endl;
}