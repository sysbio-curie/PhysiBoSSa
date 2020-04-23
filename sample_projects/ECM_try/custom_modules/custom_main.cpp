#include "custom_main.h"

using namespace PhysiCell;


std::string conc_names[2] = {"oxygen", "ecm"};

void writeDensity( int index, double t )
{
	std::string ecmname; 
	ecmname.resize( 1024 );
	sprintf( (char*) ecmname.c_str() , "_t%05d.txt", (int)round(t) );
	ecmname = "microutput//"+conc_names[index]+ecmname;
	std::ofstream outecm_file( ecmname );
	microenvironment.write_density( outecm_file, index );
	outecm_file.close();
}

/* Change the current value of the input coefficient, increase or decrease according to up value */
void evolve_coef( int up, double* coef, double dt )
/**{ 
	// increase exponentially
	if ( up )
	{
		if ( (*coef) < EPSILON ) 
			(*coef) = EPSILON; 	
		(*coef) = std::sqrt( (*coef) );
		(*coef) = (*coef) > 1 ? (1-EPSILON) : (*coef);
	}
	else
	{
		// decrease exponentially
		if ( (*coef) >= 1 )
			(*coef) = 1 - EPSILON;
		(*coef) *= (*coef);	
		(*coef) = (*coef) < 0 ? EPSILON : (*coef);
	}
}*/
{ 
	// if up, increase, else decrease
	if ( !up )
		dt = -dt;

	(*coef) +=  (*coef) * (1 - (*coef)) * dt/10.0 ;

	(*coef) = (*coef) > 1 ? (1-EPSILON) : (*coef);
	(*coef) = (*coef) < 0 ? (EPSILON) : (*coef);
}