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
#include "./Immune_system.h"
#include "./Epithelial_cell.h"
#include "./Macrophage.h"
#include "./TCell.h"

Immune_system immune_system;

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
	cell_defaults.functions.update_phenotype = NULL;  
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	int virion_index = microenvironment.find_density_index( "virion" ); 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	Cell_Definition* epitheliumCD = find_cell_definition("lung epithelium"); 
	Epithelial_Cell::setup_cell_definition(epitheliumCD);
	epitheliumCD->phenotype.sync_to_microenvironment( &microenvironment );

	immune_system.initialize();		
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	static int nV = microenvironment.find_density_index( "virion" ); 
	
	// create some cells near the origin
	
	Cell* pC;
	
	// hexagonal cell packing 
	Cell_Definition* pCD = find_cell_definition("lung epithelium"); 
	
	double cell_radius = pCD->phenotype.geometry.radius; 
	double spacing = 0.95 * cell_radius * 2.0; 
	
	double x_min = microenvironment.mesh.bounding_box[0] + cell_radius; 
	double x_max = microenvironment.mesh.bounding_box[3] - cell_radius; 

	double y_min = microenvironment.mesh.bounding_box[1] + cell_radius; 
	double y_max = microenvironment.mesh.bounding_box[4] - cell_radius; 
	
	double x = x_min; 
	double y = y_min; 
	
	double center_x = 0.5*( x_min + x_max ); 
	double center_y = 0.5*( y_min + y_max ); 
	
	double triangle_stagger = sqrt(3.0) * spacing * 0.5; 
	
	// find hte cell nearest to the center 
	double nearest_distance_squared = 9e99; 
	Cell* pNearestCell = NULL; 
	
	int n = 0; 
	while( y < y_max )
	{
		while( x < x_max )
		{
			pC = create_cell( get_cell_definition("lung epithelium" ) ); 
			pC->assign_position( x,y, 0.0 );
			
			double dx = x - center_x;
			double dy = y - center_y; 
			
			double temp = dx*dx + dy*dy; 
			if( temp < nearest_distance_squared )
			{
				nearest_distance_squared = temp;
				pNearestCell = pC; 
			}
			x += spacing; 
		}
		x = x_min; 
		
		n++; 
		y += triangle_stagger; 
		// in odd rows, shift 
		if( n % 2 == 1 )
		{
			x += 0.5 * spacing; 
		}
	}
	
	std::cout << "Infecting center cell with one virion ... " << std::endl; 
	static_cast<Epithelial_Cell*>(pNearestCell)->isInfectious = true;
	// pNearestCell->phenotype.intracellular->set_boolean_node_value("Presence_Virus", true); 
	// pNearestCell->phenotype.intracellular->set_boolean_node_value("Internal", true); 
	// pNearestCell->phenotype.intracellular->set_boolean_node_value("Replicate_Virus", true); 
	
	immune_system.create();
	
	return; 
}


std::vector<std::string> tissue_coloring_function( Cell* pCell )
{
	int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	int macrophage_type = get_cell_definition( "macrophage" ).type; 
	int tcell_type = get_cell_definition( "CD8 Tcell" ).type; 
	
	std::vector<std::string> output = {"white", "black", "white" , "white" };	

	if (pCell->type == lung_epithelial_type)
			return static_cast<Epithelial_Cell*>(pCell)->coloring_function();

	else if (pCell->type == macrophage_type)
			return static_cast<Macrophage*>(pCell)->coloring_function();
	
	else if (pCell->type == tcell_type)
			return static_cast<TCell*>(pCell)->coloring_function();
	
	return output; 
}
 
void update_tissue(double diffusion_dt)
{
	immune_system.immune_cell_recruitment( diffusion_dt ); 
	immune_system.keep_immune_cells_in_bounds( diffusion_dt ); 
}