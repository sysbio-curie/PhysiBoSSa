#ifndef __T_Cell__
#define __T_Cell__

#include "./Immune_cell.h"

using namespace BioFVM; 
using namespace PhysiCell;

class TCell : public Immune_cell
{
  private:
  
    void set_input_nodes();
    void from_nodes_to_cell();

    void function_phenotype(Phenotype& phenotype, double dt);
    void function_mechanics( Phenotype& phenotype, double dt );

    bool attempt_immune_cell_attachment(Cell* pTarget , double dt );
    Cell* check_neighbors_for_attachment(double dt);
    
  public:
    
    bool isBoundToMacrophage;
    bool isActive;
    bool isBoundToEpithelium;
    
    TCell() { isBoundToMacrophage=false; isActive=false; isBoundToEpithelium=false;}
    
    std::vector<std::string> coloring_function();
        
        
    static Cell* create_cell();
    static void setup_cell_definition(Cell_Definition* cd);
    
    inline static void function_phenotype( Cell* pCell, Phenotype& phenotype, double dt ) {
      static_cast<TCell*>(pCell)->function_phenotype(phenotype, dt);
    }
    
    inline static void function_mechanics( Cell* pCell, Phenotype& phenotype, double dt ) {
      static_cast<TCell*>(pCell)->function_mechanics(phenotype, dt);
    }
    
};

#endif