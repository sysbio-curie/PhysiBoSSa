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
	
	// Initialize MaBoSS Objects for a model
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->parse(this->network, configFile.c_str());

	IStateGroup::checkAndComplete(this->network);

	// Initialize the map relation between node name and index positions
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


/* Creates a NetworkState_Impl from the input vector */
NetworkState_Impl MaBoSSNetwork::create_networkstate(std::vector<bool>* input)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	NetworkState state;
	for (auto node : nodes)
	{
		state.setNodeState(node, (NodeState) (*input)[i]);
		i ++;
	}
	return state.getState();
}

/* Transfer state values to an output vector */
void MaBoSSNetwork::retrieve_networkstate_values(NetworkState_Impl input_state, std::vector<bool>* output)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	NetworkState state = (NetworkState) input_state;
	for ( auto node: nodes )
	{
		(*output)[i] = state.getNodeState( node ) ;
		i++;
	}
}

/* Reset a vector of bools to the init state of the network */
void MaBoSSNetwork::restart_node_values(std::vector<bool>* output)
{
	// Create MaBoSS objects to get an initial state of the network
	NetworkState network_state;
	RandomGeneratorFactory *randgen_factory = this->config->getRandomGeneratorFactory();
  	RandomGenerator *random_generator = randgen_factory->generateRandomGenerator(PhysiCell::UniformInt());
	
	// Create MaBoSS initial states
	this->network->initStates(network_state, random_generator);

	// Transfer network state to output vector
	int i = 0;
	std::vector<Node *> nodes = this->network->getNodes();
	(*output).resize(nodes.size());
	for (auto node : nodes)
	{
		(*output)[i] =  network_state.getNodeState( node );
		i++;
	}
	
	for (auto initial_value : initial_values) {
		(*node_values)[node_names[initial_value.first]] = PhysiCell::UniformRandom() < initial_value.second;
	}
}

/* Run a MaBoSS simulation with the input values*/
void MaBoSSNetwork::run_simulation(std::vector<bool>* node_values)
{	
	NetworkState_Impl state = this->create_networkstate(node_values);

	// Engine created for a single isolated simulation
	StochasticSimulationEngine* engine = new StochasticSimulationEngine(this->network, this->config);
	engine->setSeed(PhysiCell::UniformInt());
	state = engine->run(&state, NULL);
	delete engine;

	this->retrieve_networkstate_values(state, node_values);
}

/* Return the index of node based on node name */
int MaBoSSNetwork::get_maboss_node_index( std::string name )
{
	auto res = this->node_names.find(name);
	if ( res != this->node_names.end() )
		return res->second;
	
	// If node name is noy found, throw an exception
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