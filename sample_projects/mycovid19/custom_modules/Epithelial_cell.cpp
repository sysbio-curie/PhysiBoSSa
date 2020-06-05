#include "./Epithelial_cell.h" 

using namespace PhysiCell; 

std::string epithelium_submodel_version = "0.0.1"; 

Cell* Epithelial_Cell::create_cell() 
{
	return static_cast<Cell*>(new Epithelial_Cell);		
}

void Epithelial_Cell::setup_cell_definition(Cell_Definition* cd) 
{
	cd->functions.instantiate_cell = Epithelial_Cell::create_cell;
	cd->functions.update_phenotype = Epithelial_Cell::function_phenotype;
}

void Epithelial_Cell::set_input_nodes() 
{
	int virion_index = get_microenvironment()->find_density_index( "virion" );
    if (isInfected)
	{
		// std::cout << "We are infected (" << &pCell << ")" << std::endl;
		phenotype.intracellular->set_boolean_node_value("Replicate_Virus", true);
		
	} else if (nearest_density_vector()[virion_index] > custom_data["virion_detection_threshold"]) {
		// std::cout << "We are feeling the virus (" << &pCell << ")" << std::endl;
		phenotype.intracellular->set_boolean_node_value("Presence_Virus", true);
		isFeeling = true;
		// pCell->nearest_density_vector()[virion_index] -= 0.1;
	} 
	
	if (isAttachedToTCell) {
		phenotype.intracellular->set_boolean_node_value("TCellBound", true);
	}
}	

void Epithelial_Cell::from_nodes_to_cell() 
{
	if (phenotype.intracellular->get_boolean_node_value("BoundReceptor")) {
		isInContact = true;	

	}
	
	if (phenotype.intracellular->get_boolean_node_value("Replicate_Virus")) {
		isInfected = true;	
	}
	
	if (phenotype.intracellular->get_boolean_node_value("Export_Virus")) 
	{
		isInfectious = true;	

		// std::cout << "We are infectious (" << &pCell << ")" << std::endl;
		int virion_index = get_microenvironment()->find_density_index( "virion" );
		// pCell->nearest_density_vector()[virion_index] += 10;
		phenotype.secretion.net_export_rates[virion_index] = custom_data["virion_export_rate"];

	}
	
	if (phenotype.intracellular->get_boolean_node_value("Death") || phenotype.intracellular->get_boolean_node_value("DeathByTCell")) {
		static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
		start_death(apoptosis_index);

		int virion_index = get_microenvironment()->find_density_index( "virion" );
		phenotype.secretion.net_export_rates[virion_index] = 0;

	}
}

void Epithelial_Cell::function_phenotype( Phenotype& phenotype, double dt ) 
{  
	if (phenotype.intracellular->need_update()) {
		
		set_input_nodes();
		phenotype.intracellular->update();
		from_nodes_to_cell();
	}
	
	if( phenotype.death.dead == true )
	{
		functions.update_phenotype = NULL; 
		remove_all_adhesions();
	}
	
}


std::vector<std::string> Epithelial_Cell::coloring_function(  )
{
	std::vector<std::string> output( 4, "black" ); 

	if (!phenotype.death.dead)
	{
		char color [1024]; 
		int virion_index = get_microenvironment()->find_density_index( "virion" );
		// int virion_index = get_microenvironment()->find_density_index( "pro-inflammatory cytokine" );
		unsigned int gradient = (unsigned int) (1300*(nearest_density_vector()[virion_index]));
		gradient = gradient < 0 ? 0 : gradient;
		gradient = gradient > 130 ? 130 : gradient;
		
		sprintf( color, "rgb(%u,%u,%u)" , 0, 255-gradient,0 );
		// if (isFeeling)
		// 	sprintf( color, "rgb(%u,%u,%u)" , 0,125,0 );
		
		// std::cout << "Density : " << nearest_density_vector()[virion_index] << ", color = " << color << std::endl;
		if (isInContact)
			sprintf( color, "rgb(%u,%u,%u)" , 255,125,0 );
		
		if (isInfected)
			sprintf( color, "rgb(%u,%u,%u)" , 255,0,0 );
		
		if (isInfectious)
			sprintf( color, "rgb(%u,%u,%u)" , 255,125,255 );
		
		if (isAttachedToTCell)
			sprintf( color, "rgb(%u,%u,%u)" , 0,0,255 );
		
		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}
