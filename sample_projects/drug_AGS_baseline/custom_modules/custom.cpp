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
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	
	// set molecular properties 
	// int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 
	// cell_defaults.phenotype.molecular.fraction_released_at_death[tnf_substrate_index] = 0.0;
	int akti_substrate_index = microenvironment.find_density_index( "AKTi" ); 
	cell_defaults.phenotype.molecular.fraction_released_at_death[akti_substrate_index] = 0.0;
	int taki_substrate_index = microenvironment.find_density_index( "TAKi" ); 
	cell_defaults.phenotype.molecular.fraction_released_at_death[taki_substrate_index] = 0.0;
	// int bcati_substrate_index = microenvironment.find_density_index( "BCATi" ); 
	// cell_defaults.phenotype.molecular.fraction_released_at_death[bcati_substrate_index] = 0.0;


	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 

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

void update_custom_variables( Cell* pCell )
{
	// static int tnf_index = microenvironment.find_density_index( "tnf" ); 
	// static int index_tnf_concentration = pCell->custom_data.find_variable_index("tnf_concentration");
	// static int index_tnf_node = pCell->custom_data.find_variable_index("tnf_node");
	// static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
	// pCell->custom_data.variables.at(index_tnf_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[tnf_index];
	// pCell->custom_data.variables.at(index_tnf_node).value = pCell->boolean_network.get_node_value("TNF");
	// pCell->custom_data.variables.at(index_fadd_node).value = pCell->boolean_network.get_node_value("FADD");
	
	static int akti_index = microenvironment.find_density_index( "AKTi" ); 
	static int index_akti_concentration = pCell->custom_data.find_variable_index("akti_concentration");
	static int index_akti_node = pCell->custom_data.find_variable_index("akti_node");
	pCell->custom_data.variables.at(index_akti_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[akti_index];
	pCell->custom_data.variables.at(index_akti_node).value = pCell->boolean_network.get_node_value("anti_AKT");

	static int taki_index = microenvironment.find_density_index( "TAKi" ); 
	static int index_taki_concentration = pCell->custom_data.find_variable_index("taki_concentration");
	static int index_taki_node = pCell->custom_data.find_variable_index("taki_node");
	pCell->custom_data.variables.at(index_taki_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[taki_index];
	pCell->custom_data.variables.at(index_taki_node).value = pCell->boolean_network.get_node_value("anti_TAK1");

	// TNF:
	// static int nTNF_external = microenvironment.find_density_index( "tnf" );
    // double dR_EB = phenotype.molecular.internalized_total_substrates[nTNF_external];
	// static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_TNFR" ); 
		
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	// cell_defaults.custom_data[ "o2_intern" ] = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	// double o2_intern = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	// static int index_o2_intern = pCell->custom_data.find_variable_index("o2_intern");
	// pCell->custom_data.variables.at(index_o2_intern).value = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	// double o2a = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	// double o2b = phenotype.molecular.internalized_total_substrates[o2_index];
	// static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	// double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 
	double o2_in = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	double o2_ex = (pCell->nearest_density_vector())[o2_index]; // PhysiCell_constants::oxygen_index]; 
	// std::cout << "O2_Env: " << o2_ex << ". O2_Int_a: " << o2_in << "\n" << std::endl;
	// std::cout << "O2_Env: " << pO2 << ". O2_Int_a: " << o2a <<". O2_Int_b: " << o2b << "\n" << std::endl;

	// static int index_o2_intern = pCell->custom_data.find_variable_index("o2_intern");
	// static int index_o2_extern = pCell->custom_data.find_variable_index("o2_extern");
	// // o2 internal:
	// pCell->custom_data.variables.at(index_o2_intern).value = pCell->phenotype.molecular.internalized_total_substrates[o2_index];
	// // o2 external:
	// pCell->custom_data.variables.at(index_o2_extern).value = (pCell->nearest_density_vector())[o2_index];
	
	// static int bcati_index = microenvironment.find_density_index( "BCATi" ); 
	// static int index_bcati_concentration = pCell->custom_data.find_variable_index("bcati_concentration");
	// static int index_bcati_node = pCell->custom_data.find_variable_index("bcati_node");
	// pCell->custom_data.variables.at(index_bcati_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[bcati_index];
	// pCell->custom_data.variables.at(index_bcati_node).value = pCell->boolean_network.get_node_value("anti_betacatenin");

	// cell_defaults.custom_data[ "dirich_value" ] = microenvironment.dirichlet_value_vectors[0][0]
	// pCell->phenotype.molecular.internalized_total_substrates[o2_index];

}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	BooleanNetwork ags_network;
	double maboss_time_step = parameters.doubles("maboss_time_step");
	ags_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = create_cell(get_cell_definition("default")); 
		pC->assign_position( x, y, z );
		// pC->set_total_volume(sphere_volume_from_radius(radius));
		
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		pC->boolean_network = ags_network;
		pC->boolean_network.restart_nodes();
		static int index_next_physibossa_run = pC->custom_data.find_variable_index("next_physibossa_run");
		pC->custom_data.variables.at(index_next_physibossa_run).value = pC->boolean_network.get_time_to_update();
		update_custom_variables(pC);
	}

	return; 
}

// custom cell phenotype function to run PhysiBoSSa when is needed
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if ( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model) 
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) << "\n" << std::endl;
	// }
	// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0)  << "\n" << std::endl;
	// }

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	static int index_next_physibossa_run = pCell->custom_data.find_variable_index("next_physibossa_run");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (PhysiCell_globals.current_time >= pCell->custom_data.variables.at(index_next_physibossa_run).value)
	{
		set_input_nodes(pCell);

		pCell->boolean_network.run_maboss();
		// Get noisy step size
		double next_run_in = pCell->boolean_network.get_time_to_update();
		pCell->custom_data.variables.at(index_next_physibossa_run).value = PhysiCell_globals.current_time + next_run_in;
		
		update_custom_variables(pCell);

		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with live coloring 
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 
	return output; 
}

void set_input_nodes(Cell* pCell) {
	// static int tnf_index = microenvironment.find_density_index( "tnf" ); 
	// static double tnf_threshold = parameters.doubles("tnf_threshold");
	static int akti_index = microenvironment.find_density_index( "AKTi" );
	static double akti_threshold = parameters.doubles("akti_threshold");
	static int taki_index = microenvironment.find_density_index( "TAKi" );
	static double taki_threshold = parameters.doubles("taki_threshold");
	// static int bcati_index = microenvironment.find_density_index( "BCATi" );
	// static double bcati_threshold = parameters.doubles("bcati_threshold");

	// if (tnf_index != -1)
	// {
	// 	double tnf_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[tnf_index];
	// 	if (tnf_cell_concentration >= tnf_threshold)
	// 		pCell->boolean_network.set_node_value("TNF", 1);
	// 	else
	// 	{
	// 		pCell->boolean_network.set_node_value("TNF", 0);
	// 	}
	// }
	if (akti_index != -1)
	{
		double akti_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[akti_index];
		if (akti_cell_concentration >= akti_threshold)
			pCell->boolean_network.set_node_value("anti_AKT", 1);
		else
		{
			pCell->boolean_network.set_node_value("anti_AKT", 0);
		}
	}
	
	if (taki_index != -1)
	{
		double taki_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[taki_index];
		if (taki_cell_concentration >= taki_threshold)
			pCell->boolean_network.set_node_value("anti_TAK1", 1);
		else
		{
			pCell->boolean_network.set_node_value("anti_TAK1", 0);
		}		
	}
	// 	if (bcati_index != -1)
	// {
	// 	double bcati_cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[bcati_index];
	// 	if (bcati_cell_concentration >= bcati_threshold)
	// 		pCell->boolean_network.set_node_value("anti_betacatenin", 1);
	// 	else
	// 	{
	// 		pCell->boolean_network.set_node_value("anti_betacatenin", 0);
	// 	}
	// }


}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	// For AGS model

	std::string prosurvival_basename = "Prosurvival_b";
	std::string antisurvival_basename = "Antisurvival_b";
	double prosurvival_value = 0.0;
	double antisurvival_value = 0.0;

	for(int i=1; i<=3; i++)
	{
		bn_index = pCell->boolean_network.get_node_index( prosurvival_basename + std::to_string(i) );
		if ( (*nodes)[bn_index] > 0)
		{
			prosurvival_value += (*nodes)[bn_index];
		}
		bn_index = pCell->boolean_network.get_node_index( antisurvival_basename + std::to_string(i) );
		if ( (*nodes)[bn_index] > 0)
		{
			antisurvival_value += (*nodes)[bn_index];
		}
	}

	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	static int apoptosis_index; 
	double multiplier = 1.0;
	double multiplier_pro = 1.0;
	double multiplier_anti = 1.0;

	// live model 
			
	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );

		// std::cout << "Pre " << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;
		// std::cout << "Pre " << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) <<  ". " << pCell->phenotype.death.rates[apoptosis_index] << ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;

		// live rate: 0.00072
		// apoptosis.transition_rate( 0, 1) = 1.0 / (8.6 * 60.0);  = 0.001938
		
		multiplier = ( prosurvival_value + 1 ) / ( antisurvival_value + 1 ) ; //[0.25, 0.33, 0.5, 0.58, 0.66, 0.75, 1, 1.5, 2, 3, 4]
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		pCell->phenotype.death.rates[apoptosis_index] = ( 1 / multiplier ) * pCell->phenotype.death.rates[apoptosis_index];

		// multiplier = pow(10, prosurvival_value) / pow(10, antisurvival_value);
		// multiplier = ( prosurvival_value + 1 ) / ( antisurvival_value + 1) ; //[0.25, 0.33, 0.5, 0.58, 0.66, 0.75, 1, 1.5, 2, 3, 4]
		// multiplier = log2 (( prosurvival_value + 1 ) / ( antisurvival_value + 1)) ; //[- 2, 0, 2]
		// if ( prosurvival_value != 0 && antisurvival_value != 0 )
		// {
			// std::cout << "Pre " << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;
			// phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
			// std::cout << "Post " << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << ". " << multiplier << "\n" << std::endl;
		// }
		
		// multiplier_pro = pow(10, prosurvival_value);
		// multiplier_anti = pow(10, antisurvival_value);
		// phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier_pro *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		// pCell->phenotype.death.rates[apoptosis_index] = multiplier_anti * pCell->phenotype.death.rates[apoptosis_index];
		
		// multiplier_pro = 1 * prosurvival_value;
		// multiplier_anti = 0.1 * antisurvival_value;
		// if ( prosurvival_value != 0 )
		// {
		// 	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier_pro *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		// }
		// if ( antisurvival_value != 0 )
		// {
		// 	pCell->phenotype.death.rates[apoptosis_index] = multiplier_anti * pCell->phenotype.death.rates[apoptosis_index];
		// }
		
		// if ( prosurvival_value != 0 )
		// {
		// 	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = prosurvival_value *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		// }
		// if ( antisurvival_value != 0 )
		// {
		// 	pCell->phenotype.death.rates[apoptosis_index] = antisurvival_value * pCell->phenotype.death.rates[apoptosis_index];
		// }
		
		// std::cout << "Post " << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) <<  ". " << pCell->phenotype.death.rates[apoptosis_index] << ". " << prosurvival_value << ", " << antisurvival_value << ". " << multiplier_pro << ", " << multiplier_anti << "\n" << std::endl;
	}

	if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	{
		std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) <<  ". " << prosurvival_value << ", " << antisurvival_value << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;

	}

	// if ( prosurvival_value > 0 )
	// 	do_proliferation( pCell, phenotype, dt );
	// else if ( antisurvival_value > 0 )
	// {
	// 	int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
	// 	pCell->start_death(apoptosis_model_index);
	// 	return;
	// }


	// if ( pCell->boolean_network.get_node_value( "Apoptosis" ) )
	// {
	// 	int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
	// 	pCell->start_death(apoptosis_model_index);
	// 	return;
	// }

	// if ( pCell->boolean_network.get_node_value( "NonACD" ) )
	// {
	// 	int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
	// 	pCell->start_death(necrosis_model_index);
	// 	return;
	// }

	// if ( pCell->boolean_network.get_node_value( "Survival" ) )
	// {
	// 	do_proliferation( pCell, phenotype, dt );
	// }
	// int tnf_substrate_index = microenvironment.find_density_index( "tnf" );
	// static double tnf_secretion = parameters.doubles("tnf_secretion_rate");
	// double tnf_secretion_rate = 0;
	// // produce some TNF
	// if ( pCell->boolean_network.get_node_value( "NFkB" ) )
	// {
	// 	tnf_secretion_rate = (tnf_secretion / microenvironment.voxels(pCell->get_current_voxel_index()).volume);

	// }
	// pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = tnf_secretion_rate;
	pCell->set_internal_uptake_constants(dt);
}

// ***********************************************************
// * NOTE: Funtion replicated from PhysiBoSS, but not used   *
// *       as we use a live cycle model instead a Ki67 model *
// ***********************************************************
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
		//if cells in basic_Ki67_cycle_model (code 1)
	// TODO adaptar rate conforme la resta Prosurv - Antisurv
		// If cells is in G0 (quiescent) switch to pre-mitotic phase
		// if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative ) {
		// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		

	// if cells in live_cells_cycle_model (code 5)
	// TODO adaptar rate conforme la resta Prosurv - Antisurv
		// if( phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
		// if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::live) {
			// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt);
		// }
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
