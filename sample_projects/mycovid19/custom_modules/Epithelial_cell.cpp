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
	cd->functions.update_phenotype = Epithelial_Cell::phenotype;
}

void Epithelial_Cell::set_input_nodes(Cell* pCell) 
{
	int virion_index = pCell->get_microenvironment()->find_density_index( "virion" );
    if (static_cast<Epithelial_Cell*>(pCell)->isInfected)
	{
		std::cout << "We are infected (" << &pCell << ")" << std::endl;
		pCell->phenotype.intracellular->set_boolean_node_value("Presence_Virus", true);
		
	} else if (pCell->nearest_density_vector()[virion_index] > 0.1) {
		std::cout << "We are feeling the virus (" << &pCell << ")" << std::endl;
		pCell->phenotype.intracellular->set_boolean_node_value("Presence_Virus", true);
		
		pCell->nearest_density_vector()[virion_index] -= 0.1;
	}
}	

void Epithelial_Cell::from_nodes_to_cell(Cell* pCell) 
{
	if (pCell->phenotype.intracellular->get_boolean_node_value("BoundReceptor")) {
		static_cast<Epithelial_Cell*>(pCell)->isInContact = true;	

	}
	
	if (pCell->phenotype.intracellular->get_boolean_node_value("Replicate_Virus")) {
		static_cast<Epithelial_Cell*>(pCell)->isInfected = true;	
	}
	
	if (pCell->phenotype.intracellular->get_boolean_node_value("Export_Virus")) 
	{
		static_cast<Epithelial_Cell*>(pCell)->isInfectious = true;	

		std::cout << "We are infectious (" << &pCell << ")" << std::endl;
		//Export virus to microenvironment, at a specific rate
		
		// int virus_external = ; 
		// phenotype.secretion.net_export_rates[virus_external] = 1;
		int virion_index = pCell->get_microenvironment()->find_density_index( "virion" );

		pCell->nearest_density_vector()[virion_index] = 1;

	}
	
	if (pCell->phenotype.intracellular->get_boolean_node_value("Death")) {
		std::cout << "We are dying (" << &pCell << ")" << std::endl;
		static int apoptosis_index = pCell->phenotype.death.find_death_model_index( "apoptosis" ); 
		pCell->start_death(apoptosis_index);
	}
}

void Epithelial_Cell::phenotype( Cell* pCell, Phenotype& phenotype, double dt ) 
{  
	if (pCell->phenotype.intracellular->need_update()) {
		
		set_input_nodes(pCell);
		pCell->phenotype.intracellular->update();
		from_nodes_to_cell(pCell);
	}
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
	}
	
}

