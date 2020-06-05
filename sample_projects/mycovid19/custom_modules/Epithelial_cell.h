#ifndef __Epithelial_Cell__
#define __Epithelial_Cell__

#include "./Attachable_cell.h"

using namespace BioFVM; 
using namespace PhysiCell;

class Epithelial_Cell : public Attachable_cell 
{
  private:    
  public:
  
    bool isFeeling;
    bool isInContact;
    bool isInfected;
    bool isInfectious;
    bool isAttachedToTCell;
    
    Epithelial_Cell() { isInfected = false; isInContact = false; isInfectious = false; isAttachedToTCell=false;}
    
    static Cell* create_cell();
    static void setup_cell_definition(Cell_Definition* cd);
        
    void set_input_nodes();
    void from_nodes_to_cell();
    void function_phenotype(Phenotype& phenotype, double dt );
    std::vector<std::string> coloring_function();

    static void function_phenotype( Cell* pCell, Phenotype& phenotype, double dt ) {
      static_cast<Epithelial_Cell*>(pCell)->function_phenotype(phenotype, dt);
    }
};

#endif