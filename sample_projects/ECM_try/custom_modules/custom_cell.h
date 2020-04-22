#ifndef __Custom_cell_h__
#define __Custom_cell_h__

/* custom class for simulation that use cells that interact with extra cellular matrix
*/

#include "../core/PhysiCell_cell.h" 
#include "../core/PhysiCell_constants.h"
using namespace PhysiCell;

class Custom_cell : public Cell {

private:

protected:

public:
	inline bool passive() { return type == PhysiCell_constants::PASSIVE_TYPE; };
	
    /** \brief Calculate repulsion and adhesion between agent and ecm at given voxel index
	 *
	 * @param index_ecm index of the ECM density in the microenv vector of densities
	 * @param index_voxel index of the current ECM voxel  */
	void add_ecm_interaction( int index_ecm, int index_voxel );

};

#endif