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

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	cell_defaults.name = "tumor cell";
	cell_defaults.type = 0; 
	
	// set default cell cycle model
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );

	// set default_cell_functions; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;

	// set the rate terms in the default phenotype 
	// first find index for a few key variables. 
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 

	// set cycle duration, apoptotic duration and rate and initially no necrosis 
	int live_index = cell_defaults.phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
	cell_defaults.phenotype.cycle.data.transition_rate(live_index, live_index) = parameters.doubles("live_phase_duration");

	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 

	cell_defaults.phenotype.secretion.secretion_rates[tnf_substrate_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[tnf_substrate_index] = parameters.doubles("tnf_uptake_rate"); 
	cell_defaults.phenotype.secretion.saturation_densities[tnf_substrate_index] = 1; 
	
	cell_defaults.phenotype.molecular.fraction_released_at_death[tnf_substrate_index] = 0.0;

	// add custom data here, if any
	cell_defaults.custom_data.add_variable("next_physibossa_run", "dimensionless", 12.0);
	cell_defaults.custom_data.add_variable("tnf_concentration", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("tnf_node", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("fadd_node", "dimensionless", 0);

	return; 
}


void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int o2_index = microenvironment.find_density_index( "oxygen" );
	double o2 = pCell->nearest_density_vector()[o2_index];

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (PhysiCell_globals.current_time >= pCell->custom_data["next_physibossa_run"])
	{
		set_input_nodes(pCell);

		pCell->boolean_network.run_maboss();
		// Get noisy step size
		double next_run_in = pCell->boolean_network.get_time_to_update();
		pCell->custom_data["next_physibossa_run"] = PhysiCell_globals.current_time + next_run_in;
		
		update_custom_variables(pCell);

		from_nodes_to_cell(pCell, phenotype, dt);
	}
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

void update_custom_variables( Cell* pCell )
{
	std::vector<bool> * nodes = pCell->boolean_network.get_nodes();

	int tnf_maboss_index = pCell->boolean_network.get_node_index("TNF");
	int fadd_maboss_index = pCell->boolean_network.get_node_index("FADD");
	static int tnf_index = microenvironment.find_density_index( "tnf" ); 
	static double tnf_threshold = parameters.doubles("tnf_threshold");

	pCell->custom_data["tnf_concentration"] = pCell->phenotype.molecular.internalized_total_substrates[tnf_index];
	pCell->custom_data["tnf_node"] = (*nodes)[tnf_maboss_index];
	pCell->custom_data["fadd_node"] = (*nodes)[fadd_maboss_index];
}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	BooleanNetwork tnf_network;
	double maboss_time_step = parameters.doubles("maboss_time_step");
	tnf_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = create_cell(); 
		pC->assign_position( x, y, z );
		// pC->set_total_volume(sphere_volume_from_radius(radius));
		
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		pC->boolean_network = tnf_network;
		pC->boolean_network.restart_nodes();
		pC->custom_data["next_physibossa_run"] = pC->boolean_network.get_time_to_update();
		update_custom_variables(pC);
	}

	return; 
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with ki67 coloring 
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 
	return output; 
}


void set_input_nodes(Cell* pCell) {
	std::vector<bool> * nodes = pCell->boolean_network.get_nodes();
	
	int tnf_maboss_index = pCell->boolean_network.get_node_index("TNF");
	static int tnf_index = microenvironment.find_density_index( "tnf" ); 
	static double tnf_threshold = parameters.doubles("tnf_threshold");

	
	if (tnf_maboss_index != -1 && tnf_index != -1)
	{
		double tnf_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[tnf_index];
		if (tnf_cell_concentration >= tnf_threshold)
			(*nodes)[tnf_maboss_index] = 1;
		else
		{
			double rate = (tnf_threshold - tnf_cell_concentration) / tnf_threshold;
			(*nodes)[tnf_maboss_index] = uniform_random() > rate; 
		}
		
	}}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	bn_index = pCell->boolean_network.get_node_index( "Apoptosis" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	bn_index = pCell->boolean_network.get_node_index( "NonACD" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
		pCell->start_death(necrosis_model_index);
		return;
	}

	bn_index = pCell->boolean_network.get_node_index( "Survival" );
	if ( bn_index != -1 && (*nodes)[bn_index])
	{
		do_proliferation( pCell, phenotype, dt );
	}


	// For model with TNF production
	bn_index = pCell->boolean_network.get_node_index( "NFkB" );
	if ( bn_index != -1 )
	{
		int tnf_substrate_index = microenvironment.find_density_index( "tnf" );
		static double tnf_secretion = parameters.doubles("tnf_secretion_rate");

		double tnf_secretion_rate = 0;
		// produce some TNF
		if ( (*nodes)[bn_index] )
		{
			tnf_secretion_rate = (tnf_secretion / microenvironment.voxels(pCell->get_current_voxel_index()).volume);

		}
		pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = tnf_secretion_rate;
		pCell->set_internal_uptake_constants(dt);
	}
}

/* Go to proliferative if needed */
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
	// If cells is in G0 (quiescent) switch to pre-mitotic phase
	if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative )
		pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt);
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
