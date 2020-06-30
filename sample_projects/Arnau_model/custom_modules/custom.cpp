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

// std::string ecm_file;
std::vector<bool>* nodes;

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

	cell_defaults.functions.cycle_model.phase_link(1,2).arrest_function = Custom_cell::wait_for_nucleus_growth;
	cell_defaults.functions.cycle_model.phase_link(0,1).arrest_function = wait_for_cell_growth;

	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	cell_defaults.phenotype.death.models[apoptosis_model_index]->phase_link(0,1).arrest_function = Custom_cell::waiting_to_remove; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	// add custom data here, if any
	cell_defaults.custom_data.add_variable("ecm_contact", "dimensionless", 0.0); //for paraview visualization
	cell_defaults.custom_data.add_variable("pintegrin", "dimensionless", 0.5); //for paraview visualization
	cell_defaults.custom_data.add_variable("padhesion", "dimensionless", 0.5); //for paraview visualization
	cell_defaults.custom_data.add_variable("cell_contact", "dimensionless", 0.0); //for paraview visualization
	cell_defaults.custom_data.add_variable("node", "dimensionless", 0.0 ); //for paraview visualization
	build_ecm_shape();

	//Setting the custom_create_cell pointer to our create_custom_cell
	cell_defaults.functions.instantiate_cell = Custom_cell::create_custom_cell;
	// cell_defaults.functions.custom_cell_rule = Custom_cell::check_passive;
	cell_defaults.functions.update_velocity = Custom_cell::custom_update_velocity;
	cell_defaults.functions.custom_adhesion = Custom_cell::custom_adhesion_function;
	// cell_defaults.functions.add_cell_basement_membrane_interactions = Custom_cell::add_cell_basement_membrane_interactions;	
	// cell_defaults.functions.calculate_distance_to_membrane = Custom_cell::distance_to_membrane;
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
	Custom_cell* pC;
	//std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles("config_radius");
	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,tumor_radius);

	for (int i = 0; i < positions.size(); i++)
	{
		/*
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;
*/
		pC = static_cast<Custom_cell*>(create_cell(cell_defaults));
		pC->assign_position( positions[i] );
		double volume = sphere_volume_from_radius(cell_radius);
		pC->set_total_volume(volume);
		pC->phenotype.volume.target_solid_nuclear = cell_defaults.phenotype.volume.target_solid_nuclear;
		pC->phenotype.volume.target_solid_cytoplasmic = cell_defaults.phenotype.volume.target_solid_cytoplasmic;
		pC->phenotype.volume.rupture_volume = cell_defaults.phenotype.volume.rupture_volume;
		
		//pC->phenotype.cycle.data.current_phase_index = phase+1;
		//pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		//if ((phase+1) == 1)
			//pC->phenotype.cycle.pCycle_Model->phases[1].entry_function(pC, pC->phenotype, 0);

		pC->custom_data["ecm_contact"] = pC->ecm_contact;
		color_node(pC);
	}
	std::cout << "tissue created" << std::endl;

	return; 
}

std::vector<std::string> ECM_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "black" );
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	double ecm_value = pCustomCell->ecm_contact;
	int color = (int) round( ecm_value * 255.0 / (pCell->phenotype.geometry.radius)  );
	if(color > 255){
		color = 255;
	}
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255-color );
	output[0].assign( szTempString );
	return output;

}

std::vector<std::string> pMotility_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	char szTempString [128];
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	int color = (int) round(pCustomCell->pmotility * 255);
	sprintf( szTempString , "rgb(%u,0,%u)", color, 255-color );
	output[0].assign( szTempString );
	return output;
}

std::vector<std::string> node_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	if ( !pCell->phenotype.intracellular->get_boolean_node_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(0,0,255)";
		output[2] = "rgb(0,0,125)";
		
	}
	else{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	return output;
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{
	 //Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	 int color_number = parameters.ints("color_function");

	 if (color_number == 0)
	 	return ECM_coloring_function(pCell);
	 if (color_number == 1)
	 	return pMotility_coloring_function(pCell);
	 else 
	 	return node_coloring_function( pCell );
}

void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (pCell->phenotype.intracellular->need_update())
	{
		
		set_input_nodes(pCustomCell);

		pCell->phenotype.intracellular->update();
		
		from_nodes_to_cell(pCustomCell, phenotype, dt);
		color_node(pCustomCell);
	}
	pCustomCell->custom_data["ecm_contact"] = pCustomCell->ecm_contact;
	pCustomCell->custom_data["cell_contact"] = pCustomCell->cell_contact;
	
}

void set_input_nodes(Custom_cell* pCell) 
{	
	if ( pCell->phenotype.intracellular->has_node( "Oxy" ) ){
		pCell->phenotype.intracellular->set_boolean_node_value("Oxy", !pCell->necrotic_oxygen());
	}

	// 	nodes[ind] = ( !pCell->necrotic_oxygen() );

	//enough_to_node( pCell, "TGFbR", "tgfb" );

	if ( pCell->phenotype.intracellular->has_node( "Neighbours" ) ){
		pCell->phenotype.intracellular->set_boolean_node_value("Neighbours", pCell->has_neighbor(0));	
	}
	
	if ( pCell->phenotype.intracellular->has_node( "Nei2" ) ){
		pCell->phenotype.intracellular->set_boolean_node_value("Nei2", pCell->has_neighbor(1));
	}
	
	// If has enough contact with ecm or not
	if ( pCell->phenotype.intracellular->has_node( "ECM_sensing" ) )
		pCell->phenotype.intracellular->set_boolean_node_value("ECM_sensing", touch_ECM(pCell));
	
	// If has enough contact with TGFbeta or not
	if ( pCell->phenotype.intracellular->has_node( "TGFbeta" ) )
		pCell->phenotype.intracellular->set_boolean_node_value("TGFbeta", touch_TGFbeta(pCell));

	// If nucleus is deformed, probability of damage
	// Change to increase proba with deformation ? + put as parameter
	if ( pCell->phenotype.intracellular->has_node( "DNAdamage" ) )
		pCell->phenotype.intracellular->set_boolean_node_value("DNAdamage", 
			( pCell->nucleus_deform > 0.5 ) ? (2*PhysiCell::UniformRandom() < pCell->nucleus_deform) : 0
		);
		
	/// example
}

void from_nodes_to_cell(Custom_cell* pCell, Phenotype& phenotype, double dt)
{
		
	if ( pCell->phenotype.intracellular->has_node( "Apoptosis" ) 
		&& pCell->phenotype.intracellular->get_boolean_node_value( "Apoptosis" ) 
	)
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}
	
	if ( pCell->phenotype.intracellular->has_node( "Autophagy" ) 
		&& pCell->phenotype.intracellular->get_boolean_node_value( "Autophagy" ) 
	)
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	if ( pCell->phenotype.intracellular->has_node( "Hypoxia" ) 
		&& pCell->phenotype.intracellular->get_boolean_node_value( "Hypoxia" ) 
	)
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}


	if ( pCell->phenotype.intracellular->has_node( "Migration" ))
		pCell->set_oxygen_motility(pCell->phenotype.intracellular->get_boolean_node_value("Migration"));
		
		// pCell->evolve_motility_coef( pCell->phenotype.intracellular->get_boolean_node_value( "Migration" ), dt );

	 if ( pCell->phenotype.intracellular->has_node( "Cell_growth" ) && pCell->phenotype.intracellular->get_boolean_node_value("Cell_growth") ){
	 	//do_proliferation( pCell, phenotype, dt );
	 }

	/*if ( pCell->phenotype.intracellular->has_node("CCA") 
		&& pCell->phenotype.intracellular->get_boolean_node_value("CCA") )
	{
		pCell->freezing(1);
	}	
	*/

	/*if ( pCell->phenotype.intracellular->has_node( "Polarization" ) )
		pCell->evolve_polarity_coef( 
			pCell->phenotype.intracellular->get_boolean_node_value( "Polarization" ), dt 
		);
	*/

	if ( pCell->phenotype.intracellular->has_node( "Cell_cell" ) )
		pCell->evolve_cellcell_coef( 
			pCell->phenotype.intracellular->get_boolean_node_value( "Cell_cell" ), dt 
		);

	if ( pCell->phenotype.intracellular->has_node( "Matrix_adhesion" ) )
		pCell->evolve_integrin_coef( 
			pCell->phenotype.intracellular->get_boolean_node_value( "Matrix_adhesion" ), dt 
		);

	if ( pCell->phenotype.intracellular->has_node("Matrix_modif") )
	 	pCell->set_mmp( pCell->phenotype.intracellular->get_boolean_node_value("Matrix_modif") );


	if ( pCell->phenotype.intracellular->has_node("EMTreg") )
	{
		pCell->set_mmp( pCell->phenotype.intracellular->get_boolean_node_value("EMTreg") );
	}

	pCell->freezing( 0 );

	if ( pCell->phenotype.intracellular->has_node( "Quiescence" ) 
		&& pCell->phenotype.intracellular->get_boolean_node_value( "Quiescence" )
	){
		pCell->freezing(1);
		cell_defaults.functions.volume_update_function = static_volume_function;
	}
	
	if ( pCell->phenotype.intracellular->has_node( "Quiescence" ) 
		&& !pCell->phenotype.intracellular->get_boolean_node_value( "Quiescence" )
	){
		pCell->freezing(1);
		cell_defaults.functions.volume_update_function = standard_volume_update_function;
	}
	if ( pCell->phenotype.intracellular->has_node( "Cell_freeze" ) ){
		pCell->freezer(3 * pCell->phenotype.intracellular->get_boolean_node_value( "Cell_freeze" ));
	}
}




void build_ecm_shape() {
 
 // Here we design a spherical shell of ecm
 std::vector<double> center(3, 0);
 double inner_radius = parameters.doubles("config_radius");
 //double outer_radius = 150;
 
 for (auto voxel : microenvironment.mesh.voxels) {
 // Compute norm to center
 double t_norm = norm(voxel.center);
 // If norm is in [inner_radius, outer_radius], then we add it
 if (t_norm > inner_radius ) {
 microenvironment.density_vector(voxel.mesh_index)[microenvironment.find_density_index("ecm")] = 0.3; 
 microenvironment.density_vector(voxel.mesh_index)[microenvironment.find_density_index("TGFbeta")] = 0.3;
 }
 }
 
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

/* Go to proliferative if needed */
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
	// If cells is in G0 (quiescent) switch to pre-mitotic phase
	if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative )
		pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt);
}

double sphere_volume_from_radius(double rad)
{
	double PI4_3 = 4.0 / 3.0 * M_PI;
return PI4_3 * rad * rad * rad;
}

bool touch_ECM(Custom_cell* pCell)
{ 
	return pCell->contact_ecm() > parameters.doubles("contact_cell_ECM_threshold"); 
}

bool touch_TGFbeta(Custom_cell* pCell){
	return pCell->contact_TGFbeta() > parameters.doubles("ECM_TGFbeta_ratio");
}

void enough_to_node( Custom_cell* pCell, std::string nody, std::string field )
{
	if ( pCell->phenotype.intracellular->has_node(nody) )
	{
		int felt = pCell->feel_enough(field, *pCell);
		if ( felt != -1 )
			pCell->phenotype.intracellular->set_boolean_node_value(nody, felt);
	}
}

void color_node(Custom_cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data["node"] = pCell->phenotype.intracellular->get_boolean_node_value(node_name);
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}

bool wait_for_cell_growth(Cell* pCell, Phenotype& phenotype, double dt){

	return !pCell->phenotype.intracellular->get_boolean_node_value("Cell_growth");

}