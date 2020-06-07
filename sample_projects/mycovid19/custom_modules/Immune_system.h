#ifndef __Immune_system__
#define __Immune_system__

#include "../core/PhysiCell.h"

class Immune_system
{
  private:    
    void create_infiltrating_Tcell(void);
    void keep_immune_cells_off_edge( void );

  public:
    void initialize();
    void create();
    void immune_cell_recruitment( double dt );
    void keep_immune_cells_in_bounds( double dt );
  
};



#endif