#include <vector>
#include "../core/PhysiCell_cell.h"
#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_mesh.h"
#include "../BioFVM/BioFVM_microenvironment.h"

using namespace PhysiCell;

class Custom_cell;

class Custom_cell_container : public Cell_Container
{
private:
	std::vector<Custom_cell*> cells_ready_to_divide; // the index of agents ready to divide
	std::vector<Custom_cell*> cells_ready_to_die;
    bool initialzed = false;
public:
    void update_all_cells(double t, double phenotype_dt, double mechanics_dt, double diffusion_dt ); 
    void custom_flag_cell_for_division( Custom_cell* pCell ); 
	void custom_flag_cell_for_removal( Custom_cell* pCell ); 
};