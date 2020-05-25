#ifndef __Custom_cell_h__
#define __Custom_cell_h__

/* custom class for simulation that use cells that interact with extra cellular matrix
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include "../core/PhysiCell_cell.h" 
#include "../core/PhysiCell_constants.h"
#include "custom_main.h"
using namespace PhysiCell;

class Sphere : public Cell {

private:
	
protected:
	
public:
	inline bool passive() { return type == PhysiCell_constants::PASSIVE_TYPE; };
	
    Sphere();

};

#endif