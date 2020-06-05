#ifndef __Immune_cell__
#define __Immune_cell__

#include "./Attachable_cell.h"

using namespace BioFVM; 
using namespace PhysiCell;

class Immune_cell : public Attachable_cell
{
  private:    
  public:

    void keep_in_bounds(double Xmin, double Xrange, double Ymin, double Yrange, double Zmin, double Zrange);
    
};



#endif