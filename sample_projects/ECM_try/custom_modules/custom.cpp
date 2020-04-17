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
#include "../BioFVM/BioFVM.h"  
#include "../addons/PhysiBoSSa/src/boolean_network.h"

using namespace BioFVM;

// declare cell definitions here 

// std::string ecm_file;
std::vector<bool> nodes;

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

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// add custom data here, if any
	cell_defaults.custom_data.add_variable("next_physibossa_run", "dimensionless", 12.0);
	std::cout << microenvironment.find_density_index("ecm");
	
	load_ecm_file();

	// set the rate terms in the default phenotype 
	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int ecm_substrate_index = microenvironment.find_density_index("ecm");


	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 

	//cell_defaults.phenotype.secretion.uptake_rates[ecm_substrate_index] = 10; 
	//cell_defaults.phenotype.secretion.secretion_rates[ecm_substrate_index] = 0; 
	//cell_defaults.phenotype.secretion.saturation_densities[ecm_substrate_index] = 38; 
	
	microenvironment.diffusion_coefficients[ecm_substrate_index] = 1e-85;
	microenvironment.decay_rates[ecm_substrate_index] = 0;
	microenvironment.list_indexes(0.5);
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
	MaBoSSNetwork* maboss;
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	BooleanNetwork ecm_network;
	ecm_network.initialize_boolean_network(bnd_file, cfg_file, 12);

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
		//std::cout<<(*all_cells)[i]->position<<std::endl;
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		pC->boolean_network = ecm_network;
		pC->boolean_network.restart_nodes();
		pC->custom_data["next_physibossa_run"] = pC->boolean_network.get_time_to_update();
		//std::cout<< pC->position.size() << std::endl;
		//std::cout<< pC->position << std::endl;
	}
	std::cout<<(*all_cells)[25]<<std::endl;
	std::cout << "tissue created" << std::endl;

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with live dead coloring 
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 
	return output; 
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
		
		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

void set_input_nodes(Cell* pCell) {
int ind;
	nodes = *(pCell->boolean_network.get_nodes());
	// Oxygen input node O2; Oxygen or Oxy
	ind = pCell->boolean_network.get_node_index( "Oxygen" );
	if ( ind < 0 )
		ind = pCell->boolean_network.get_node_index( "Oxy" );
	if ( ind < 0 )
		ind = pCell->boolean_network.get_node_index( "O2" );
	if ( ind >= 0 )
		nodes[ind] = ( !pCell->necrotic_oxygen() );
	

	ind = pCell->boolean_network.get_node_index( "Neighbours" );
	if ( ind >= 0 )
		nodes[ind] = ( pCell->has_neighbor(0) );
	
	ind = pCell->boolean_network.get_node_index( "Nei2" );
	if ( ind >= 0 )
		nodes[ind] = ( pCell->has_neighbor(1) );

	// If has enough contact with ecm or not
	ind = pCell->boolean_network.get_node_index( "ECM_sensing" );
	if ( ind >= 0 )
		nodes[ind] = ( parameters.ints("contact_cell_ECM_threshold") );
	// If has enough contact with ecm or not
	ind = pCell->boolean_network.get_node_index( "ECM" );
	if ( ind >= 0 )
		nodes[ind] = ( parameters.ints("contact_cell_ECM_threshold") );
	// If has enough contact with ecm or not
	ind = pCell->boolean_network.get_node_index( "ECMicroenv" );
	if ( ind >= 0 )
		nodes[ind] = ( parameters.ints("contact_cell_ECM_threshold") );
	
	// If nucleus is deformed, probability of damage
	// Change to increase proba with deformation ? + put as parameter
	ind = pCell->boolean_network.get_node_index( "DNAdamage" );
	//std::cout << mycell->nucleus_deformation() << std::endl;
	if ( ind >= 0 )
		nodes[ind] = ( pCell->nucleus_deform > 0.5 ) ? (2*PhysiCell::UniformRandom() < pCell->nucleus_deform) : 0;
	/// example
}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* point_to_nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	bn_index = pCell->boolean_network.get_node_index( "Apoptosis" );
	if ( bn_index != -1 && (*point_to_nodes)[bn_index] )
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	bn_index = pCell->boolean_network.get_node_index("Matrix_modif");
	if ( bn_index != -1 && (*point_to_nodes)[bn_index] )
	{
		pCell->set_mmp( (*point_to_nodes)[bn_index] );
		return;
	}

	bn_index = pCell->boolean_network.get_node_index("EMT");
	if ( bn_index != -1 && (*point_to_nodes)[bn_index] )
	{
		pCell->set_mmp( (*point_to_nodes)[bn_index] );
		return;
	}

	/// example
}


/* Load ecm density values from given file */
void load_ecm_file()
{
	// strip( &ecm_file );
	std::cout << "Loading ECM file " << parameters.strings("init_ecm_filename") << std::endl;
	std::ifstream infile;
	infile.open( parameters.strings("init_ecm_filename") );
	std::ofstream outfile ("verify_output.txt");
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
		
		outfile << "voxel_index: " << voxel_index << "   " << "ecm_index: " << microenvironment.find_density_index("ecm") << "   " << "amount_of_density: " << microenvironment.density_vector(voxel_index)[microenvironment.find_density_index("ecm")] << std::endl;
		
	}
	outfile.close();
	infile.close();
	std::cout << "File loaded !" << std::endl;

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
