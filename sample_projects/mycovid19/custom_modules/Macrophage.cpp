#include "./Macrophage.h" 

using namespace PhysiCell; 

Cell* Macrophage::create_cell() 
{
	return static_cast<Cell*>(new Macrophage);		
}

void Macrophage::setup_cell_definition(Cell_Definition* cd) 
{
	cd->functions.instantiate_cell = Macrophage::create_cell;
	cd->functions.update_phenotype = Macrophage::function_phenotype;
}

void Macrophage::set_input_nodes() 
{
	int virion_index = get_microenvironment()->find_density_index( "virion" );
	hasDetectedVirus = nearest_density_vector()[virion_index] > user_parameters->doubles("macrophage_virion_detection_threshold");
	phenotype.intracellular->set_boolean_node_value("Presence_Virus", hasDetectedVirus);
}	

void Macrophage::from_nodes_to_cell() 
{
	int virion_index = get_microenvironment()->find_density_index( "virion" );
	int cytokines_index = get_microenvironment()->find_density_index( "pro-inflammatory cytokine" );
	
	isActive = phenotype.intracellular->get_boolean_node_value("Active");
	if (isActive){
		phenotype.secretion.uptake_rates[virion_index] = user_parameters->doubles("macrophage_eating_rate");
		phenotype.secretion.secretion_rates[cytokines_index] = user_parameters->doubles("macrophage_cytokin_release_rate");
	} else {
		phenotype.secretion.uptake_rates[virion_index] = 0;
		phenotype.secretion.secretion_rates[cytokines_index] = 0;
	}
}

void Macrophage::function_phenotype(Phenotype& phenotype, double dt ) 
{  
	if (this->phenotype.intracellular->need_update()) {
		
		this->set_input_nodes();
		this->phenotype.intracellular->update();
		this->from_nodes_to_cell();
	}
	
	if( this->phenotype.death.dead == true )
	{
		this->functions.update_phenotype = NULL; 
		remove_all_adhesions();
	}
	
}

std::vector<std::string> Macrophage::coloring_function()
{
	std::vector<std::string> output( 4, "black" ); 

	if( phenotype.death.dead == false )
	{
		char color [1024]; 
		
		int virion_index = get_microenvironment()->find_density_index( "virion" );
		unsigned int gradient = (unsigned int) ((130/user_parameters->doubles("macrophage_virion_detection_threshold"))*(nearest_density_vector()[virion_index]));
		gradient = gradient < 0 ? 0 : gradient;
		gradient = gradient > 130 ? 130 : gradient;
	
		sprintf( color, "rgb(%u,%u,%u)" ,0,gradient, 255 );

		if (isActive)
			sprintf( color, "rgb(%u,%u,%u)", 0, 255, 255);
			
		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}
