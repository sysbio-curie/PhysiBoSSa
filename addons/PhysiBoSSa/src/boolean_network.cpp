#include "boolean_network.h"

void BooleanNetwork::initialize_boolean_network(std::string bnd_file, std::string cfg_file, double time_step){
	this->maboss.init_maboss(bnd_file, cfg_file);
	this->maboss.set_update_time_step(time_step);
}

void BooleanNetwork::restart_nodes() 
{
	this->maboss.restart_node_values(&(this->nodes));
	this->set_time_to_update();
}

/* Update MaboSS network states */
void BooleanNetwork::run_maboss()
{
	this->maboss.run_simulation(&this->nodes);
	this->set_time_to_update();
}