#include "MaBoSSNetwork.h"

/* Default constructor */
MaBoSSNetwork::MaBoSSNetwork( std::string networkFile, std::string configFile )
{
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->setSeedPseudoRandom( UniformInt() );
	this->config->parse(this->network, configFile.c_str());

	//This can be modified
	this->update_time_step = 10;

	IStateGroup::checkAndComplete(network);

	std::vector<Node *> nodes = this->network->getNodes();
	int i = 0;
	for (auto node : nodes)
	{
		this->node_names[ node->getLabel() ] = i;
		i++;	
	}

	StochasticSimulationEngine* engine = new StochasticSimulationEngine(this->network, this->config);
	int seed = this->config->getSeedPseudoRandom();
	engine->setSeed(seed);

	(*this->state) = engine->run(NULL, NULL);
	delete engine;
}

/* Default destructor */
MaBoSSNetwork::~MaBoSSNetwork()
{
	delete this->network;
	this->network = NULL;
	delete this->config;
	this->config = NULL;
}

/* Run the current network */
void MaBoSSNetwork::run(std::vector<bool>* nodes_val)
{
	StochasticSimulationEngine* engine = new StochasticSimulationEngine(this->network, this->config);
	int seed = this->config->getSeedPseudoRandom();
	engine->setSeed(seed);

	(*this->state) = engine->run(this->state, NULL);

	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		(*nodes_val)[i] = ((NetworkState) (*this->state)).getNodeState( node ) ;
		i++;
	}
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes()
{
	int i = 0;
	std::vector<Node*> nodes = network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << ((NetworkState) (*this->state)).getNodeState( node ) << "; ";
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