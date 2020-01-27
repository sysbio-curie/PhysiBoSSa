#include "MaBoSSCellCycleNetwork.h"

CellCycleNetwork::CellCycleNetwork(std::string bnd_file, std::string cfg_file)
{
	MaBoSSNetwork* maboss = new MaBoSSNetwork(bnd_file, cfg_file); 
	this->set_maboss(maboss);
	this->time_to_update = ( 1 + 0.5*UniformRandom11() ) * maboss->get_update_time_step();
}

CellCycleNetwork::~CellCycleNetwork()
{
	delete this->maboss;
	this->maboss = NULL;
}

/* Initialization: set network */
void CellCycleNetwork::set_maboss( MaBoSSNetwork* maboss )
{
	this->maboss = maboss;
	this->set_time_to_update();

	// initialize all nodes to 0
	this->nodes.resize( maboss->nodes_lenght() );
}

/* random update time, to asynchronize it between all cells */
void CellCycleNetwork::set_time_to_update()
{
	this->time_to_update = (1 + 0.5*UniformRandom11()) * this->maboss->get_update_time_step();
}

/* Update MaboSS network states */
void CellCycleNetwork::run_maboss()
{
	this->maboss->run(&this->nodes);
	this->set_time_to_update();
}