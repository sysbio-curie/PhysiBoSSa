#include "MaBoSSNetwork.h"

/* Default constructor */
MaBoSSNetwork::MaBoSSNetwork( std::string networkFile, std::string configFile )
{
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->setSeedPseudoRandom( UniformInt() );
	this->config->parse(this->network, configFile.c_str());

	this->state = NULL;

	//This can be modified
	this->update_time_step = 10;

	IStateGroup::checkAndComplete(network);

	int i = 0;
	std::vector<Node *> nodes = this->network->getNodes();
	for (auto node : nodes)
	{
		this->node_names[ node->getLabel() ] = i;
		i++;	
	}

	this->run();
}

/* Default destructor */
MaBoSSNetwork::~MaBoSSNetwork()
{
	delete this->network;
	this->network = NULL;
	delete this->config;
	this->config = NULL;
}

void MaBoSSNetwork::load_state(std::vector<bool>* input)
{
	int i = 0;
	std::vector<Node*> nodes = network->getNodes();
	for (auto node : nodes)
	{
		((NetworkState) *(this->state)).setNodeState(node, (NodeState) (*input)[i]);
		i ++;
	}
}

void MaBoSSNetwork::recover_state(std::vector<bool>* output)
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		(*output)[i] = ((NetworkState) *(this->state)).getNodeState( node ) ;
		i++;
	}
}

/* Run the current network */
void MaBoSSNetwork::run()
{
	StochasticSimulationEngine* engine = new StochasticSimulationEngine(this->network, this->config);
	int seed = this->config->getSeedPseudoRandom();
	engine->setSeed(seed);
	this->next_state = engine->run(this->state, NULL);
	this->state = &next_state;
	delete engine;
}

/* Run the current network */
void MaBoSSNetwork::run(std::vector<bool>* node_values)
{
	this->load_state(node_values);
	this->run();
	this->recover_state(node_values);
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes()
{
	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << ((NetworkState) *(this->state)).getNodeState( node ) << "; ";
		i++;
	}
	std::cout << std::endl;
}

/* Return the index of node based on node's name */
int MaBoSSNetwork::get_node_index( std::string name )
{
	auto res = this->node_names.find(name);
	if ( res != this->node_names.end() )
		return res->second;
	return -1;	
}