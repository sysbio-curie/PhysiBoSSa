#ifndef __Attachable_cell__
#define __Attachable_cell__

#include "../core/PhysiCell.h"
#include "../core/PhysiCell_phenotype.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

class Attachable_cell : public PhysiCell::Cell 
{
  private:    
  public:
    
    void attach( Cell* pCell );
    void detach( Cell* pCell );
    void remove_all_adhesions();
  
};



#endif