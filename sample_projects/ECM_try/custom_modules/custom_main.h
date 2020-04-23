#include <string>
#include <iostream>
#include <vector>
#include "../core/PhysiCell.h"


/** \brief List of voxels indexes to write to output files */

void writeDensity( int index, double t );
/* Change the current value of the input coefficient, increase or decrease according to up value */
void evolve_coef( int up, double* coef, double dt );