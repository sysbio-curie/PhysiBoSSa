#include "./Macrophage.h" 

using namespace PhysiCell; 

std::string macrophage_submodel_version = "0.0.1"; 

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

	if (nearest_density_vector()[virion_index] > custom_data["virion_detection_threshold"]) {
		phenotype.intracellular->set_boolean_node_value("Presence_Virus", true);
		hasDetectedVirus = true;
	}
}	

void Macrophage::from_nodes_to_cell() 
{

	if (phenotype.intracellular->get_boolean_node_value("Eating")){
		int virion_index = get_microenvironment()->find_density_index( "virion" );
		phenotype.secretion.uptake_rates[virion_index] = custom_data["macrophage_eating_rate"];
	}
	
	if (phenotype.intracellular->get_boolean_node_value("Release_Cytokines")) {
		int cytokines_index = get_microenvironment()->find_density_index( "pro-inflammatory cytokine" );
		phenotype.secretion.secretion_rates[cytokines_index] = custom_data["macrophage_cytokin_release_rate"];
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
		
		sprintf( color, "rgb(%u,%u,%u)" , 0,0, 255 );
		
		if (hasDetectedVirus)
			sprintf( color, "rgb(%u,%u,%u)" , 125,125,255 );
		
		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}
