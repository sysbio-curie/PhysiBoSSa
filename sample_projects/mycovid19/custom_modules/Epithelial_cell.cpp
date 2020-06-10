#include "./Epithelial_cell.h" 

using namespace PhysiCell; 

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
	
	// If we are marked as infectious, but the boolean model is not
	// This in the case used by the initialization, to start as infectious
    if (isInfectious & !phenotype.intracellular->get_boolean_node_value("Infectious"))
	{
		phenotype.intracellular->set_boolean_node_value("Presence_Virus", true);
		phenotype.intracellular->set_boolean_node_value("Infected", true);
		phenotype.intracellular->set_boolean_node_value("Infectious", true);
		
	} else {
		bool presence_virus = nearest_density_vector()[virion_index] > user_parameters->doubles("virion_detection_threshold");
		phenotype.intracellular->set_boolean_node_value("Presence_Virus", presence_virus);
	}
	
	if (isAttachedToTCell) {
		phenotype.intracellular->set_boolean_node_value("TCellBound", true);
	}
}	

void Epithelial_Cell::from_nodes_to_cell() 
{
	isInfected = phenotype.intracellular->get_boolean_node_value("Infected");
	isInfectious = phenotype.intracellular->get_boolean_node_value("Infectious");
	
	if (isInfectious)
	{
		int virion_index = get_microenvironment()->find_density_index( "virion" );
		phenotype.secretion.net_export_rates[virion_index] = user_parameters->doubles("virion_export_rate");
	}
	
	if (phenotype.intracellular->get_boolean_node_value("Death") || phenotype.intracellular->get_boolean_node_value("DeathByTCell")) {
		remove_all_adhesions();

		int virion_index = get_microenvironment()->find_density_index( "virion" );
		phenotype.secretion.net_export_rates[virion_index] = 0;

		static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
		start_death(apoptosis_index);
	}
	
	if (phenotype.intracellular->get_boolean_node_value("CureByTCell")) {
		remove_all_adhesions();

		int virion_index = get_microenvironment()->find_density_index( "virion" );		
		phenotype.secretion.net_export_rates[virion_index] = 0;
		isCured = true;
		isInfected = false;
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
		unsigned int gradient = (unsigned int) ((255/user_parameters->doubles("virion_detection_threshold"))*(nearest_density_vector()[virion_index]));
		gradient = gradient < 0 ? 0 : gradient;
		gradient = gradient > 255 ? 255 : gradient;
		
		sprintf( color, "rgb(%u,%u,%u)" , gradient, 255,0 );
		
		if (isInfected)
			sprintf( color, "rgb(%u,%u,%u)" , 255,125,0 );
		
		if (isInfectious)
			sprintf( color, "rgb(%u,%u,%u)" , 255,0,0 );
		
		if (isCured)
			sprintf( color, "rgb(%u,%u,%u)" , 0,125, 0 );
		
		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}
