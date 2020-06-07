#ifndef __Epithelial_Cell__
#define __Epithelial_Cell__

#include "./Attachable_cell.h"

using namespace BioFVM; 
using namespace PhysiCell;

class Epithelial_Cell : public Attachable_cell 
{
  private:    
    
    void set_input_nodes();
    void from_nodes_to_cell();
    void function_phenotype(Phenotype& phenotype, double dt );
    
  public:
  
    bool isInfected;
    bool isInfectious;
    bool isAttachedToTCell;
    bool isCured;
    
    Epithelial_Cell() { isInfected = false; isInfectious = false; isAttachedToTCell=false; isCured=false;}
    
    static Cell* create_cell();
    static void setup_cell_definition(Cell_Definition* cd);
        
    std::vector<std::string> coloring_function();

    static void function_phenotype( Cell* pCell, Phenotype& phenotype, double dt ) {
      static_cast<Epithelial_Cell*>(pCell)->function_phenotype(phenotype, dt);
    }
};

#endif