#ifndef __Custom_cell_container_h__
#define __Custom_cell_container_h__

#include <vector>
#include "../core/PhysiCell_cell.h"
#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_mesh.h"
#include "../BioFVM/BioFVM_microenvironment.h"
#include <string>

using namespace PhysiCell;

class Custom_cell;

class Custom_cell_container : public Cell_Container
{
protected:
    std::vector<Custom_cell*> *custom_cells_ready_to_divide; // the index of agents ready to divide
	std::vector<Custom_cell*> *custom_cells_ready_to_die;
    double membrane_length;
    std::string membrane_shape;

public:
    Custom_cell_container();
    void update_all_cells(double t );
    void update_all_cells(double t, double phenotype_dt, double mechanics_dt, double diffusion_dt ); 
    double get_membrane_length();
    std::string get_membrane_shape();
    /** \brief Set the size of the BM */
	inline void set_membrane_length( double l )
	{ membrane_length = l; } ;
    /** \brief Set the shape of the BM */
	inline void set_membrane_shape( std::string shape )
	{ membrane_shape = shape; strip( &membrane_shape ); } ;

};

extern std::vector<Custom_cell*> *all_custom_cells;
Custom_cell_container* create_custom_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size );
#endif