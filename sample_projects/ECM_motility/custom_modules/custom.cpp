/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include "./custom.h"
#include "../BioFVM/BioFVM.h"  

using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> * nodes;

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling; 
	cell_defaults.functions.instantiate_cell = Custom_cell::create_custom_cell;

	cell_defaults.functions.custom_cell_rule = Custom_cell::check_passive;
	cell_defaults.functions.update_velocity = Custom_cell::custom_update_velocity;
	cell_defaults.functions.custom_adhesion = Custom_cell::custom_adhesion_function;

	cell_defaults.functions.cycle_model.phase_link(1,2).arrest_function = Custom_cell::wait_for_nucleus_growth;

	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	cell_defaults.phenotype.death.models[apoptosis_model_index]->phase_link(0,1).arrest_function = Custom_cell::waiting_to_remove; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	// add custom data here, if any

	load_ecm_file();

	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	Custom_cell* pC;
	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	
	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = static_cast<Custom_cell*>(create_cell(cell_defaults));
		pC->assign_position( x, y, z );
		pC->set_initial_volume(cell_defaults, radius);
		pC->phenotype.cycle.data.current_phase_index = phase+1;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		if ((phase+1) == 1)
			pC->phenotype.cycle.pCycle_Model->phases[1].entry_function(pC, pC->phenotype, 0);
		
	}
	std::cout << "tissue created" << std::endl;

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	std::vector< std::string > output( 4 , "black" );
	double ecm_value = pCustomCell->ecm_contact;
	int color = (int) round( (ecm_value * 255)/ (pCell->phenotype.geometry.radius) );
	if (color > 255)
	{
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255-color );
	output[0].assign( szTempString );

	//std::vector< std::string > output = false_cell_coloring_live_dead(pCell);
	return output; 
}

void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (pCell->phenotype.intracellular->need_update())
	{
		set_input_nodes(pCell);

		pCell->phenotype.intracellular->update();
				
		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

void set_input_nodes(Cell* pCell) 
{
	if ( pCell->phenotype.intracellular->has_node("ECMicroenv") ){
		Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
		pCell->phenotype.intracellular->set_boolean_node_value("ECMicroenv", touch_ECM(pCustomCell));
	}
	// If nucleus is deformed, probability of damage
	// Change to increase proba with deformation ? + put as parameter
	if ( pCell->phenotype.intracellular->has_node("ECMicroenv") ){
	
		Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);

		pCell->phenotype.intracellular->set_boolean_node_value("DNAdamage", 
			(pCustomCell->nucleus_deform > 0.5 ) ? 
				(2*PhysiCell::UniformRandom() < pCustomCell->nucleus_deform) : 0
		);
	}
}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	if (pCell->phenotype.intracellular->has_node( "Apoptosis" )
	 	&& pCell->phenotype.intracellular->get_boolean_node_value("Apoptosis") 
	)
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	if ( pCell->phenotype.intracellular->has_node( "Migration" ) )
		static_cast<Custom_cell*>(pCell)->evolve_motility_coef( 
			pCell->phenotype.intracellular->get_boolean_node_value("Migration"), dt 
		);
	
	if (pCell->phenotype.intracellular->has_node("CCA")
		&& pCell->phenotype.intracellular->get_boolean_node_value("CCA")
	)
		static_cast<Custom_cell*>(pCell)->freezing(1);
	
	if (pCell->phenotype.intracellular->has_node("EMT")
		&& pCell->phenotype.intracellular->get_boolean_node_value("EMT")
	)
		static_cast<Custom_cell*>(pCell)->set_mmp(
			pCell->phenotype.intracellular->get_boolean_node_value("EMT") 
		);
}


/* Load ecm density values from given file */
void load_ecm_file()
{
	// strip( &ecm_file );
	std::cout << "Loading ECM file " << parameters.strings("init_ecm_filename") << std::endl;
	std::ifstream infile;
	infile.open( parameters.strings("init_ecm_filename") );
	std::string array[4];
	int i = 0;
	std::string line;
	//skip first line: title
	getline( infile, line, '\n' ); 
	while ( getline( infile, line, '\n') )
	{
		std::stringstream ss;
		ss.str( line );
		i = 0;
		while ( getline( ss, array[i], ';') )
		{
			i++;
		}
		double x = std::stod(array[0]);
		double y = std::stod(array[1]);
		double z = std::stod(array[2]);
		double amount = std::stod(array[3]);
		std::vector<double> pos(3);
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		int voxel_index = microenvironment.nearest_voxel_index( pos );
		microenvironment.density_vector(voxel_index)[microenvironment.find_density_index("ecm")] += amount; 		
	}
	infile.close();
}


std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header) 
{ 
	// File pointer 
	std::fstream fin; 
	std::vector<init_record> result;

	// Open an existing file 
	fin.open(filename, std::ios::in); 

	// Read the Data from the file 
	// as String Vector 
	std::vector<std::string> row; 
	std::string line, word;

	if(header)
		getline(fin, line);

	do 
	{
		row.clear(); 

		// read an entire row and 
		// store it in a string variable 'line' 
		getline(fin, line);

		// used for breaking words 
		std::stringstream s(line); 

		// read every column data of a row and 
		// store it in a string variable, 'word' 
		while (getline(s, word, delimiter)) { 

			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}

bool touch_ECM(Custom_cell* pCell)
{ 
	return pCell->contact_ecm() > parameters.doubles("contact_cell_ECM_threshold"); 
}