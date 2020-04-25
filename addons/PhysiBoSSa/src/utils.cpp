#include "utils.h"
#include "../../../core/PhysiCell_utilities.h"

double UniformRandom11()
{
	double res = PhysiCell::UniformRandom();
	return ( 2.0 * ( res - 0.5 ) ); 
}
