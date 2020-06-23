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

#include "./custom.h"

#include "./tnf_response_dynamics.h"

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 

		
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	std::cout << cell_defaults.name<< std::endl;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment ); 

	// custom parameter
	/*
	Parameter<double> paramD;

	paramD = parameters.doubles["TNFR_binding_rate"];
	cell_defaults.custom_data[ "TNFR binding rate" ] = paramD.value;

	paramD = parameters.doubles["TNFR_endocytosis_rate"];
	cell_defaults.custom_data[ "TNFR endocytosis rate" ] = paramD.value;

	paramD = parameters.doubles["TNFR_recycling_rate"];
	cell_defaults.custom_data[ "TNFR recycling rate" ] = paramD.value;

	paramD = parameters.doubles["TNFR_receptors_per_cell"]; 
	cell_defaults.custom_data[ "unbound external TNFR" ] = paramD.value;

	paramD = parameters.doubles["TFN_net_production_rate"]; 
	cell_defaults.custom_data[ "TFN net production rate" ] = paramD.value;

	paramD = parameters.doubles["TNFR_activation_threshold"]; 
	cell_defaults.custom_data[ "TNFR activation threshold" ] = paramD.value;
	
	*/

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();
	std::cout << cell_defaults.name << std::endl;
	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	tnf_dynamics_model_setup();
	submodel_registry.display( std::cout ); 
		
	// set molecular properties 
	int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 
	cell_defaults.phenotype.molecular.fraction_released_at_death[tnf_substrate_index] = 0.0;

	float tnf_uptake_rates = cell_defaults.custom_data[ "TNFR binding rate" ] * parameters.doubles["TNFR_receptors_per_cell"].value;
	cell_defaults.phenotype.secretion.uptake_rates[tnf_substrate_index] = 0;
	cell_defaults.custom_data[ "unbound_external_TNFR" ] = cell_defaults.custom_data["TNFR_receptors_per_cell"]; 

	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	

	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{

	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	double maboss_time_step = parameters.doubles("maboss_time_step");

	BooleanNetwork tnf_network;
	tnf_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);

	for (int i = 0; i < cells.size(); i++)
	{
		Cell* pC;
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = create_cell(get_cell_definition("default")); 
		pC->assign_position( x, y, z );
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		pC->boolean_network = tnf_network;
		pC->boolean_network.restart_nodes();

		static int index_next_physibossa_run = pC->custom_data.find_variable_index("next_physibossa_run");
		pC->custom_data[index_next_physibossa_run] = pC->boolean_network.get_time_to_update();
		update_monitor_variables(pC);
	}

	return; 
}

// custom cell phenotype function to run PhysiBoSSa when is needed

void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	
	static int index_next_physibossa_run = pCell->custom_data.find_variable_index("next_physibossa_run");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (PhysiCell_globals.current_time >= pCell->custom_data[index_next_physibossa_run])
	{
		
		tnf_dynamics_model(pCell, phenotype, dt );
		update_boolean_model_input(pCell, phenotype, dt );
		pCell->boolean_network.run_maboss();
		update_cell_state_model_based(pCell, phenotype, dt);

		// Get noisy step size
		double next_run_in = pCell->boolean_network.get_time_to_update();
		pCell->custom_data[index_next_physibossa_run] = PhysiCell_globals.current_time + next_run_in;
	}
	update_monitor_variables(pCell);
}


void update_monitor_variables(Cell* pCell )
{
    //static int index_tnf_external = microenvironment.find_density_index( "tnf" ); 
	//static int tnf_internal = pCell->custom_data.find_variable_index( "tnf" );
    
	static int index_tnf_node = pCell->custom_data.find_variable_index("tnf_node");
	static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
	static int index_nfkb_node = pCell->custom_data.find_variable_index("nfkb_node");
	
	pCell->custom_data[index_nfkb_node] = pCell->boolean_network.get_node_value( "NFkB" ) ;
	pCell->custom_data[index_tnf_node] = pCell->boolean_network.get_node_value("TNF");
	pCell->custom_data[index_fadd_node] = pCell->boolean_network.get_node_value("FADD");

	// pCell->custom_data[index_tnf_external]= microenvironment.nearest_density_vector;
	// pCell->custom_data[tnf_internal];	
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with live coloring 
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 

	// dead cells 
	if( pCell->phenotype.death.dead == false )
	{
		static int nR_EB = pCell->custom_data.find_variable_index( "bound external TNFR" );  
		float activation_threshold = pCell->custom_data.find_variable_index( "TNFR activation threshold" );

		int bounded_tnf = (int) round( (pCell->custom_data[nR_EB] / activation_threshold) * 255.0 ); 
		if (bounded_tnf > 0)
		{
			char szTempString [128];
			sprintf( szTempString , "rgb(%u,%u,%u)", bounded_tnf, bounded_tnf, 255-bounded_tnf );
			output[0].assign( "black" );
			output[1].assign( szTempString );
			output[2].assign( "black" );
			output[3].assign( szTempString );
		}
		//output[0] = "blue"; 
		//output[2] = "darkblue"; 
		//output[1] = "red"; 
		//output[3] = "darkred";
	}
	
	

	
	
	return output;
}
// ***********************************************************
// * NOTE: Funtion to read init files created with PhysiBoSS *
// ***********************************************************
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
