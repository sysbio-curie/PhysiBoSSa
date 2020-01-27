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
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = Ki67_advanced;
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 
	// first find index for a few key variables. 
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 

	// set cycle duration, apoptotic duration and rate and initially no necrosis 
	int ki67n_index = cell_defaults.phenotype.cycle.model().find_phase_index(PhysiCell_constants::Ki67_negative);
	int ki67pre_index = cell_defaults.phenotype.cycle.model().find_phase_index(PhysiCell_constants::Ki67_positive_premitotic);
	int ki67post_index = cell_defaults.phenotype.cycle.model().find_phase_index(PhysiCell_constants::Ki67_positive_postmitotic);
	
	cell_defaults.phenotype.cycle.data.transition_rate(ki67n_index, ki67pre_index) = parameters.doubles("ki67n_phase_duration");
	cell_defaults.phenotype.cycle.data.transition_rate(ki67pre_index, ki67post_index) = parameters.doubles("ki67pre_phase_duration");
	cell_defaults.phenotype.cycle.data.transition_rate(ki67post_index, ki67n_index) = parameters.doubles("ki67post_phase_duration");

	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 

	cell_defaults.phenotype.secretion.uptake_rates[tnf_substrate_index] = parameters.doubles("tnf_uptake_rate"); 
	cell_defaults.phenotype.secretion.secretion_rates[tnf_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[tnf_substrate_index] = 1; 
	
	// add custom data here, if any
	cell_defaults.custom_data.add_variable("next_physibossa_run", "dimensionless", 10.0);
	cell_defaults.functions.custom_cell_rule = boolean_network_rule;

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
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
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
		pC->set_total_volume(sphere_volume_from_radius(radius));
		
		pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		pC->maboss_cycle_network = new CellCycleNetwork(bnd_file, cfg_file);
		pC->custom_data["next_physibossa_run"] = pC->maboss_cycle_network->get_time_to_update();
	}

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with ki67 coloring 
	std::vector<std::string> output = false_cell_coloring_Ki67(pCell); 
	return output; 
}

void boolean_network_rule(Cell* pCell, Phenotype& phenotype, double dt )
{
	if (PhysiCell_globals.current_time > pCell->custom_data["next_physibossa_run"] && !phenotype.death.dead) 
	{
		set_input_nodes(pCell);

		#pragma omp critical
		pCell->maboss_cycle_network->run_maboss();
		
		double next_run_in = pCell->maboss_cycle_network->get_time_to_update();
		pCell->custom_data["next_physibossa_run"] = PhysiCell_globals.current_time + next_run_in;
		
		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

void set_input_nodes(Cell* pCell) {
	std::vector<bool> * nodes = pCell->maboss_cycle_network->get_nodes();
	
	int tnf_maboss_index = pCell->maboss_cycle_network->get_maboss_node_index("TNF");
	static int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 
	static double tnf_threshold = parameters.doubles("tnf_threshold");

	if (tnf_maboss_index != -1 && tnf_substrate_index != -1)
		(*nodes)[tnf_maboss_index] = (*pCell->bounded)[tnf_substrate_index] > tnf_threshold;
}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->maboss_cycle_network->get_nodes();
	static double tnf_secretion = parameters.doubles("tnf_secretion_rate");

	int bn_index = pCell->maboss_cycle_network->get_maboss_node_index( "Survival" );
	if ( bn_index != -1 && (*nodes)[bn_index])
		do_proliferation( pCell, phenotype, dt );

	bn_index = pCell->maboss_cycle_network->get_maboss_node_index( "Apoptosis" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		delete pCell->maboss_cycle_network;
		pCell->maboss_cycle_network = NULL;
		pCell->functions.update_phenotype = NULL;
		return;
	}

	bn_index = pCell->maboss_cycle_network->get_maboss_node_index( "NonACD" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
		pCell->start_death(necrosis_model_index);
		delete pCell->maboss_cycle_network;
		pCell->maboss_cycle_network = NULL;
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// For model with TNF production
	bn_index = pCell->maboss_cycle_network->get_maboss_node_index( "NFkB" );
	if ( bn_index != -1 )
	{
		int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 
		// produce some TNF
		pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = 
			(*nodes)[bn_index] ? (tnf_secretion / microenvironment.voxels(pCell->get_current_voxel_index()).volume) : 0;
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