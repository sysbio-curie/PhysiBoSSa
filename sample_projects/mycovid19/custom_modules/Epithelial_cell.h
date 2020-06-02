#ifndef __Epithelial_Cell__
#define __Epithelial_Cell__

#include "../core/PhysiCell.h"
#include "../core/PhysiCell_phenotype.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

class Epithelial_Cell : public PhysiCell::Cell 
{
  private:    
  public:
  
    bool isInContact;
    bool isInfected;
    bool isInfectious;
    
    Epithelial_Cell() { isInfected = false; isInContact = false; isInfectious = false;}
    
    static Cell* create_cell();
    static void setup_cell_definition(Cell_Definition* cd);
        
    static void set_input_nodes(Cell* pCell);
    static void from_nodes_to_cell(Cell* pCell);
    
    static void phenotype( Cell* pCell, Phenotype& phenotype, double dt );
};

#endif