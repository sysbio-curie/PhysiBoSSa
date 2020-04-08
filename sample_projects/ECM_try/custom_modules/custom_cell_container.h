#ifndef __Custom_cell_container_h__
#define __Custom_cell_container_h__

#include <vector>
#include "../core/PhysiCell_cell.h"
#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_mesh.h"
#include "../BioFVM/BioFVM_microenvironment.h"

using namespace PhysiCell;

class Custom_cell;

class Custom_cell_container : public Cell_Container
{
protected:
    std::vector<Custom_cell*> *custom_cells_ready_to_divide; // the index of agents ready to divide
	std::vector<Custom_cell*> *custom_cells_ready_to_die;

public:
    Custom_cell_container();
    void update_all_cells(double t );
    void update_all_cells(double t, double phenotype_dt, double mechanics_dt, double diffusion_dt ); 

};

extern std::vector<Custom_cell*> *all_custom_cells;
#endif