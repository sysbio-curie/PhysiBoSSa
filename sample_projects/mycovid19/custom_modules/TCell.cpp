#include "./TCell.h" 
#include "./Macrophage.h"
#include "./Epithelial_cell.h"

using namespace PhysiCell; 

Cell* TCell::create_cell() 
{
	return static_cast<Cell*>(new TCell);		
}

void TCell::setup_cell_definition(Cell_Definition* cd) 
{
	cd->functions.instantiate_cell = TCell::create_cell;
	cd->functions.update_phenotype = TCell::function_phenotype;
	cd->functions.custom_cell_rule = TCell::function_mechanics;
}

void TCell::set_input_nodes() 
{	
	phenotype.intracellular->set_boolean_node_value("Contact_Macrophage", isBoundToMacrophage);
	phenotype.intracellular->set_boolean_node_value("Contact_Epithelium", isBoundToEpithelium);
}	

void TCell::from_nodes_to_cell() 
{	
	if (this->phenotype.intracellular->get_boolean_node_value("Active")) {
		isActive = true;
		remove_all_adhesions();

		// Activated TCell will hunt for infected cells which export virus. So virus is used for chemotaxis
		phenotype.motility.chemotaxis_index = get_microenvironment()->find_density_index( "virion");
		// bias direction is gradient for the indicated substrate 
		phenotype.motility.migration_bias_direction = nearest_gradient(phenotype.motility.chemotaxis_index);
		// move up or down gradient based on this direction 
		phenotype.motility.migration_bias_direction *= phenotype.motility.chemotaxis_direction; 

		// normalize 
		normalize( &( phenotype.motility.migration_bias_direction ) );
	}
}

void TCell::function_phenotype(Phenotype& phenotype, double dt ) 
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


bool TCell::attempt_immune_cell_attachment(Cell* pTarget , double dt )
{
	// Checking if the is in our targets, otherwise return		
	if (
		pTarget->type == get_cell_definition( "CD8 Tcell" ).type 
		
	||
		(pTarget->type == get_cell_definition( "macrophage" ).type 
		&& !(static_cast<Macrophage*>(pTarget)->isActive && !isActive))
		
	|| 	(pTarget->type == get_cell_definition( "lung epithelium").type
		&& !(static_cast<Epithelial_Cell*>(pTarget)->isInfected && !static_cast<Epithelial_Cell*>(pTarget)->isCured && isActive)) 
	) 
	{ return false; }

	// if the target cell is dead, give up 
	if( pTarget->phenotype.death.dead == true )
	{ return false; } 

	// if the target cell is too far away, give up 
	std::vector<double> displacement = pTarget->position - position;
	double distance_scale = norm( displacement ); 
	if( distance_scale > user_parameters->doubles("max_attachment_distance") )
	{ return false; } 

	// now, get the attachment probability 
	double attachment_probability = user_parameters->doubles("cell_attachment_rate") * dt; 

	// don't need to cap it at 1.00: if prob > 100%, 
	// then this statement always evaluates as true, 
	// just the same as capping probability at 100% 
	if( UniformRandom() <= attachment_probability )
	{
		attach(pTarget);
		
		if (pTarget->type == get_cell_definition( "macrophage" ).type) {
			isBoundToMacrophage = true;
			static_cast<Macrophage*>(pTarget)->isAttachedToTCell = true;
			
		} else if (pTarget->type == get_cell_definition( "lung epithelium" ).type) {
			isBoundToEpithelium = true;
			static_cast<Epithelial_Cell*>(pTarget)->isAttachedToTCell = true;
		}
		return true; 
	}
	
	return false; 	
}


Cell* TCell::check_neighbors_for_attachment(double dt)
{
	std::vector<Cell*> nearby = cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != this )
		{
			if( attempt_immune_cell_attachment( nearby[i] , dt ) )
			{
				return nearby[i]; 
			}
		}
		i++; 
	}
	
	return NULL; 
}



void TCell::function_mechanics(Phenotype& phenotype, double dt)
{
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// detach all attached cells 
		remove_all_adhesions(); 
		
		// Let's just fully disable now. 
		functions.custom_cell_rule = NULL; 
		functions.update_phenotype = NULL;

		return; 
	}	
	
	// if I am not adhered to a cell, turn motility on 
	if( state.neighbors.size() == 0 )
	{ phenotype.motility.is_motile = true; }
	else
	{ phenotype.motility.is_motile = false; return; }	
		
	// I'm not attached, look for cells nearby and try to attach
	
	// if this returns non-NULL, we're now attached to a cell 
	if( check_neighbors_for_attachment( dt) )
	{
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	
	isBoundToEpithelium = false;
	isBoundToMacrophage = false;
	
	return; 
}


std::vector<std::string> TCell::coloring_function()
{
	std::vector<std::string> output( 4, "black" ); 

	if( phenotype.death.dead == false )
	{
		char color [1024]; 
		
		sprintf( color, "rgb(%u,%u,%u)" ,255 ,255, 255 );
		
		if (isBoundToMacrophage)
			sprintf( color, "rgb(%u,%u,%u)" , 125,125,255 );
		
		if (isActive)
			sprintf( color, "rgb(%u,%u,%u)" , 0,125,255 );
		
		if (isBoundToEpithelium)
			sprintf( color, "rgb(%u,%u,%u)" , 255,0,0 );
		
		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
	
	return output; 
}
