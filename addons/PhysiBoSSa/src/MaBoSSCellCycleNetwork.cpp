#include "MaBoSSCellCycleNetwork.h"

void CellCycleNetwork::restart_nodes() 
{
	this->maboss.restart_node_values(&(this->nodes));
	this->set_time_to_update();
}

/* Update MaboSS network states */
void CellCycleNetwork::run_maboss()
{
	this->maboss.run_simulation(&this->nodes);
	this->set_time_to_update();
}