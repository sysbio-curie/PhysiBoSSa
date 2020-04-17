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
