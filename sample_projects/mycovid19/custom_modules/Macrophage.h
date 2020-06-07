#ifndef __Macrophage__
#define __Macrophage__

#include "./Immune_cell.h"

using namespace BioFVM; 
using namespace PhysiCell;

class Macrophage : public Immune_cell
{
  private:    
  public:
  
    bool hasDetectedVirus;
    bool isActive;
    bool isAttachedToTCell;
    Macrophage() { hasDetectedVirus = false; isActive=false; isAttachedToTCell=false;}
    
    static Cell* create_cell();
    static void setup_cell_definition(Cell_Definition* cd);
        
    void set_input_nodes();
    void from_nodes_to_cell();
    void function_phenotype(Phenotype& phenotype, double dt);
    std::vector<std::string> coloring_function();

    static void function_phenotype( Cell* pCell, Phenotype& phenotype, double dt ) {
      static_cast<Macrophage*>(pCell)->function_phenotype(phenotype, dt);
    }
    
};

#endif