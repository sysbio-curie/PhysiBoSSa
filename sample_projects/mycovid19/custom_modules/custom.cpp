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
#include "./Epithelial_cell.h"
#include "./Macrophage.h"
#include "./TCell.h"

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

	Cell_Definition* macrophagesCD = find_cell_definition("macrophage"); 
	Macrophage::setup_cell_definition(macrophagesCD);
	macrophagesCD->phenotype.sync_to_microenvironment( &microenvironment );
	
	Cell_Definition* tcellCD = find_cell_definition("CD8 Tcell"); 
	TCell::setup_cell_definition(tcellCD);
	tcellCD->phenotype.sync_to_microenvironment( &microenvironment );
			
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
	static_cast<Epithelial_Cell*>(pNearestCell)->isInfected = true; 
	
	// now place immune cells 
	for (int i=0; i < parameters.ints("number_macrophages"); i++) {
	
		pC = create_cell( get_cell_definition("macrophage" ) ); 
		pC->assign_position( 
			x_min + (x_max-x_min)*UniformRandom(),
			y_min + (y_max-y_min)*UniformRandom(),
			0.0 + parameters.doubles("immune_z_offset")
		);
		
	}	
		
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
 
void keep_immune_cells_off_edge( void )
{
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 

	static bool setup_done = false; 
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// warning hardcoded
	static double relative_edge_margin = 0; // 0.1; 
	static double relative_interior = 1 - 2 * relative_edge_margin; 
	
	if( setup_done == false )
	{
		Xmin += relative_edge_margin*Xrange; 
		Ymin += relative_edge_margin*Yrange; 
		Zmin += relative_edge_margin*Zrange;
		
		Xrange *= relative_interior;
		Yrange *= relative_interior;
		Zrange *= relative_interior;  
		setup_done = true; 
	}
	
	static int epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	
	for( int n=0 ; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if (pC->type != epithelial_type)
			static_cast<Immune_cell*>(pC)->keep_in_bounds(Xmin, Xrange, Ymin, Yrange, Zmin, Zrange);
	}
	return; 
}


void keep_immune_cells_in_bounds( double dt )
{
	static double dt_bounds = 5; 
	static double next_time = 0.0; 

	static double t_bounds = 0.0; 
	static double t_last_bounds = 0.0; 
	static double t_next_bounds = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_bounds > t_next_bounds- tolerance )
	{
		double elapsed_time = (t_bounds - t_last_bounds );
		
		keep_immune_cells_off_edge(); 
		
		t_last_bounds = t_bounds; 
		t_next_bounds = t_bounds + dt_bounds; 
	}
	t_bounds += dt; 

	return; 
}


void create_infiltrating_Tcell(void)
{
	// static Cell_Definition* pCD = find_cell_definition( "CD8 Tcell" );
	// create_infiltrating_immune_cell( pCD ); 

	// return; 
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static bool setup_done = false; 
	
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// keep cells away from the outer edge 
	
	if( setup_done == false )
	{
		Xmin += 0.1*Xrange; 
		Ymin += 0.1*Yrange; 
		Zmin = 0;
		
		Xrange *= 0.8;
		Yrange *= 0.8;
		Zrange = 0.0; 
		setup_done = true; 
	}
	
	// create some of each type of cell 
	
	Cell* pC;
	
	std::vector<double> position = {0,0,0}; 
	position[0] = Xmin + UniformRandom()*Xrange; 
	position[1] = Ymin + UniformRandom()*Yrange; 
	position[2] = parameters.doubles("immune_z_offset"); 
		
	pC = create_cell( get_cell_definition("CD8 Tcell" ) ); 
	pC->assign_position( position );
		
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	pC->phenotype.secretion.uptake_rates[proinflammatory_cytokine_index] = 
		parameters.doubles("activated_cell_cytokine_uptake_rate"); // 10;
	
	return;
}



void immune_cell_recruitment( double dt )
{
	static int proinflammatory_cytokine_index = 
		microenvironment.find_density_index("pro-inflammatory cytokine");
	
	static double dt_immune = parameters.doubles( "immune_dt" ); 
	static double t_immune = 0.0; 
	static double t_last_immune = 0.0; 
	static double t_next_immune = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_immune > t_next_immune- tolerance )
	{
		double elapsed_time = (t_immune - t_last_immune );
//		std::cout << "Immune time! " << t_immune << " (elapsed: " << elapsed_time << ") " << std::endl; 
		
		
		// CD8 T cell recruitment 
		
		static double CD8_Tcell_recruitment_rate = parameters.doubles( "CD8_Tcell_max_recruitment_rate" ); 
		static double TC_min_signal = parameters.doubles( "CD8_Tcell_recruitment_min_signal" ); 
		static double TC_sat_signal = parameters.doubles( "CD8_Tcell_recruitment_saturation_signal" ); 
		static double TC_max_minus_min = TC_sat_signal - TC_min_signal; 
		
		double total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
		double total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[proinflammatory_cytokine_index] - TC_min_signal ); 
			dRate /= TC_max_minus_min; 
			// crop to [0,1] 
			if( dRate > 1 ) 
			{ dRate = 1; } 
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate; 
		}	
		// multiply by dV and rate_max 
		total_scaled_signal = total_rate; 
		
		total_rate *= microenvironment.mesh.dV; 
		total_rate *= CD8_Tcell_recruitment_rate; 
		
		// expected number of new neutrophils 
		int number_of_new_cells = (int) round( total_rate * elapsed_time ); 
		if( number_of_new_cells )
		{
			// std::cout << "\tRecruiting " << number_of_new_cells << " CD8 T cells ... " << std::endl; 
			
//			std::cout << "\tTotal signal/dV : " << total_scaled_signal << std::endl;
//			std::cout << "\tTotal signa : " << total_scaled_signal * microenvironment.mesh.dV << std::endl; 
//			double total_volume = microenvironment.mesh.dV * microenvironment.mesh.voxels.size() ; 
//			std::cout << "\tmean signal : " << total_scaled_signal * microenvironment.mesh.dV / total_volume << std::endl; 
			for( int n = 0; n < number_of_new_cells ; n++ )
			{ create_infiltrating_Tcell(); }
		}
		
		t_last_immune = t_immune; 
		t_next_immune = t_immune + dt_immune; 
		
		// std::cout << "\t\tnext immune time: " << t_next_immune << std::endl;  
	}
	t_immune += dt; 
	return; 
}

void update_tissue(double diffusion_dt)
{
	immune_cell_recruitment( diffusion_dt ); 
	
	keep_immune_cells_in_bounds( diffusion_dt ); 
}