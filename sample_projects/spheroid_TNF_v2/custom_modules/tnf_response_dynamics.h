
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

void tnf_dynamics_model_setup();

void tnf_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt );

void update_boolean_model_input( Cell* pCell, Phenotype& phenotype, double dt );

void update_cell_state_model_based(Cell* pCell, Phenotype& phenotype, double dt);

