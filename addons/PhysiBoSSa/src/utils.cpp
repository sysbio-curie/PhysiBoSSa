#include "utils.h"
#include "../core/PhysiCell_utilities.h"

std::random_device rd;
std::mt19937 gen(rd());

int UniformInt()
{
	std::uniform_int_distribution<int> int_dis;
	return int_dis(gen);
}

double UniformRandom11()
{
	double res = PhysiCell::UniformRandom();
	return ( 2.0 * ( res - 0.5 ) ); 
}
