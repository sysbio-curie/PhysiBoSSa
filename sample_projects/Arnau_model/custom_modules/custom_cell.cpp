#include "custom_cell.h"

Custom_cell::Custom_cell() {
	freezed = 0;
    pintegrin = 0.5;
	pmotility = 0.5;
	padhesion = 0.5;
	ecmrad = sqrt(3.0) * get_microenvironment()->mesh.dx / 2.0;
	ecm_contact = 0;
	TGFbeta_contact = 0;
	cell_contact = PhysiCell::parameters.doubles("initial_cell_contact_parameter");
	nucleus_deform = 0;
}

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void Custom_cell::add_ecm_interaction( int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		displacement = position - get_container()->underlying_mesh.voxels[index_voxel].center;
		double distance = norm(displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);
		
		double dd = phenotype.geometry.radius + ecmrad;  
		double dnuc = phenotype.geometry.nuclear_radius + ecmrad;  

		double tmp_r = 0;
		// Cell overlap with ECM node, add a repulsion term
		if ( distance < dd )
		{
			// repulsion stronger if nucleii overlap, see Macklin et al. 2012, 2.3.1
			if ( distance < dnuc )
			{
				double M = 1.0;
				double c = 1.0 - dnuc/dd;
				c *= c;
				c -= M;
				tmp_r = c*distance/dnuc + M;
				nucleus_deform += (dnuc-distance);
			}
			else
			{
				tmp_r = ( 1 - distance / dd );
				tmp_r *= tmp_r;
			}
			tmp_r *= dens * PhysiCell::parameters.doubles("cell_ecm_repulsion");
		}

		// Cell adherence to ECM through integrins
		double max_interactive_distance = (PhysiCell::parameters.doubles("max_interaction_factor")*phenotype.geometry.radius) + ecmrad;
		if ( distance < max_interactive_distance ) 
		{	
			double temp_a = 1 - distance/max_interactive_distance; 
			temp_a *= temp_a; 
			/* \todo change dens with a maximal density ratio ? */
			ecm_contact += dens * (max_interactive_distance-distance);
			// temp_a *= dens * ( static_cast<Cell*>(this) )->integrinStrength();
			temp_a *= dens * integrinStrength();
			tmp_r -= temp_a;
		}
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		velocity += tmp_r * displacement;
	}
}

/** \brief (De)-Activate ECM degradation by the cell */
void Custom_cell::set_mmp( int activate )
{ 
	int ecm_index = get_microenvironment()->find_density_index("ecm");

	if (activate)
		phenotype.secretion.uptake_rates[ecm_index] = PhysiCell::parameters.doubles("ecm_degradation") * pintegrin;
		
	else
		phenotype.secretion.uptake_rates[ecm_index] = 0;
	
}

/* Return value of adhesion strength with ECM according to integrin level */
double Custom_cell::integrinStrength()
{ 
	return get_integrin_strength( pintegrin ); 
}

/* Return if cell has enough contact with other cells (compared to given threshold determined by the given level) */	
bool Custom_cell::has_neighbor(int level)
{ 
	if ( level == 0 ){
		//std::cout << contact_cell() << " - " << PhysiCell::parameters.doubles("contact_cell_cell_threshold") << std::endl;
		return contact_cell() > PhysiCell::parameters.doubles("contact_cell_cell_threshold"); 
		}
	else{
		//std::cout << contact_cell() << " - " << PhysiCell::parameters.doubles("contact_cell_cell_threshold") << std::endl;
		return contact_cell() > (2 * PhysiCell::parameters.doubles("contact_cell_cell_threshold")); 
	}
}

/* Calculate adhesion coefficient with other cell */
double Custom_cell::adhesion( Cell* other_cell )
{
    Custom_cell* custom_other_cell = static_cast<Custom_cell*>(other_cell);
	
	double adh = 0;

	if (PhysiCell::parameters.ints("choose_adhesion_function") == 0){
		//bool my_cell_cell = phenotype.intracellular->get_boolean_node_value("Cell_cell");
		//bool your_cell_cell = custom_other_cell->phenotype.intracellular->get_boolean_node_value("Cell_cell");
		bool my_single = phenotype.intracellular->get_boolean_node_value("Single");
		bool your_single = custom_other_cell->phenotype.intracellular->get_boolean_node_value("Single");
		bool my_mig = phenotype.intracellular->get_boolean_node_value("Migration");
		bool your_mig = custom_other_cell->phenotype.intracellular->get_boolean_node_value("Migration");
		if (my_single || your_single )
			return adh;
		else if ( my_mig && your_mig || !my_mig && !your_mig){
			adh = PhysiCell::parameters.doubles("homotypic_adhesion_max") * padhesion;
			return adh;
			}	
		else if(!my_mig && your_mig || my_mig && !your_mig){
			adh = 0.3 * padhesion;
			return adh;
		}
		/*
		else if(my_mig && !your_mig){
			adh = 0.3 * padhesion;
			return adh;
		}
		
		else if (!my_cell_cell || !your_cell_cell){
			adh = std::min( get_homotypic_strength(padhesion), custom_other_cell->get_homotypic_strength(padhesion) );
			return adh;
		}
		
		else if (my_cell_cell && !my_mig && your_cell_cell && !your_mig){
			adh = sqrt(get_homotypic_strength(padhesion) * custom_other_cell->get_homotypic_strength(padhesion));
			return adh;
		}
		*/
	}
	else{
		if(phenotype.intracellular->get_boolean_node_value("Single") || custom_other_cell->phenotype.intracellular->get_boolean_node_value("Single"))
			return adh;
	
		if ( this->type == other_cell->type )
			adh = std::min( get_homotypic_strength(padhesion), custom_other_cell->get_homotypic_strength(padhesion) );
		else
			adh = std::min( get_heterotypic_strength(padhesion), custom_other_cell->get_heterotypic_strength(padhesion) );
	}
	
	return adh;
}


double Custom_cell::get_adhesion()
{
	if (passive())
		return 0;
	else
		return 1;
}

void Custom_cell::set_oxygen_motility(bool active)
{	
	//phenotype.motility.is_motile = active;
		
	if (active){
		phenotype.motility.chemotaxis_index = get_microenvironment()->find_density_index( "oxygen");
		// bias direction is gradient for the indicated substrate 
		phenotype.motility.migration_bias_direction = nearest_gradient(phenotype.motility.chemotaxis_index);
		phenotype.motility.migration_bias = PhysiCell::parameters.doubles("migration_bias");
		phenotype.motility.chemotaxis_direction = 1.0;
		phenotype.motility.migration_speed = PhysiCell::parameters.doubles("migration_speed");
		phenotype.motility.persistence_time = PhysiCell::parameters.doubles("persistence");
		// move up or down gradient based on this direction 
		phenotype.motility.migration_bias_direction *= phenotype.motility.chemotaxis_direction * get_motility_amplitude(pmotility); 
	}
	else{
		//restore to default
		phenotype.motility.chemotaxis_index = cell_defaults.phenotype.motility.chemotaxis_index; 
		phenotype.motility.migration_bias_direction = cell_defaults.phenotype.motility.migration_bias_direction;
		phenotype.motility.migration_bias = cell_defaults.phenotype.motility.migration_bias;
		phenotype.motility.chemotaxis_direction = cell_defaults.phenotype.motility.chemotaxis_direction;
		phenotype.motility.migration_speed = cell_defaults.phenotype.motility.migration_speed;
		phenotype.motility.persistence_time = cell_defaults.phenotype.motility.persistence_time;
	 	}
};

/* Update the value of freezing of the cell with bitwise operation
* Do a bitwise-or comparison on freezed and input parameter:
* if freezed = 0, it will be the value of the parameter frozen
* if freezed = 1, it will be either 1 (frozen = 0) or 3 (frozen = 3) */
void Custom_cell::freezer( int frozen )
{
    freezed = freezed | frozen;
}

Cell* Custom_cell::create_custom_cell()
{
	Custom_cell* pNew; 
	pNew = new Custom_cell;		
	return static_cast<Cell*>(pNew); 
}

void Custom_cell::custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);

	pCustomCell->ecm_contact = 0;
	pCustomCell->nucleus_deform = 0;
	pCustomCell->TGFbeta_contact = 0;
	pCustomCell->cell_contact = 0;
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	for( auto neighbor : pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()] )
	{
		pCell->add_potentials( neighbor );
	}

	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 )
		pCustomCell->add_ecm_interaction( ecm_index, pCell->get_current_mechanics_voxel_index() );
		pCustomCell->add_TGFbeta_interaction(pCell->get_current_mechanics_voxel_index());

	for (auto neighbor_voxel_index : pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()])
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[neighbor_voxel_index].center, neighbor_voxel_index))
			continue;

		if ( ecm_index >= 0 )
			pCustomCell->add_ecm_interaction( ecm_index, neighbor_voxel_index );
			pCustomCell->add_TGFbeta_interaction(pCell->get_current_mechanics_voxel_index());
		for( auto other_neighbor : pCell->get_container()->agent_grid[neighbor_voxel_index] )
		{
			pCell->add_potentials(other_neighbor);
		}
	}
	
	//std::cout << pCustomCell->cell_contact << "  ";
	/*
	if((pCustomCell->ecm_contact) > (pCustomCell->cell_contact * PhysiCell::parameters.doubles("ecm_cell_contact_factor"))){
		pCustomCell->padhesion = 0;
	}
*/
	if (pCustomCell->freezed > 2){
		return ;
	}

	
		pCell->update_motility_vector(dt);
		//std::cout << phenotype.motility.motility_vector << "  ";
		pCell->velocity += phenotype.motility.motility_vector;
	
	return; 
}

double Custom_cell::custom_repulsion_function(Cell* pCell, Cell* otherCell) 
{
	if ((static_cast<Custom_cell*>(pCell))->passive())
		return 0; // Sphere are not affected by repulsion
	else
		return sqrt( pCell->phenotype.mechanics.cell_cell_repulsion_strength * otherCell->phenotype.mechanics.cell_cell_repulsion_strength ); 
}	
	
double Custom_cell::custom_adhesion_function(Cell* pCell, Cell* otherCell, double distance) 
{

	Custom_cell* custom_pCell = static_cast<Custom_cell*>(pCell);
	Custom_cell* custom_otherCell = static_cast<Custom_cell*>(otherCell);
	// hyp: junction strength is limited by weakest cell
	double adh;
	double thisadh = custom_pCell->get_adhesion();
	double otadh = custom_otherCell->get_adhesion();

	// first case, passive cell with active cell
	if ( thisadh == 0 && otadh == 1 )
	{
		custom_pCell->ecm_contact += distance;
		adh = custom_otherCell->integrinStrength();
	}
	else
	{
		// second case, active cell with passive cell
		if ( thisadh == 1 && otadh == 0 )
		{
			custom_pCell->ecm_contact += distance;
			adh = custom_pCell->integrinStrength();
		}
		else
		{
			// passive, passive
			if ( thisadh == 0 && otadh == 0 )
			{
				adh = 0;
			}
			// active, active
			else
			{
				//std::cout << distance << "  ";
				custom_pCell->cell_contact += distance;
				adh = custom_pCell->adhesion(custom_otherCell);
				custom_pCell->custom_data["adh"] = adh;
			}
		}
	}
	return adh;
}

bool Custom_cell::wait_for_nucleus_growth (Cell* cell, Phenotype& phenotype, double dt) {
	return (relative_diff( 
		phenotype.volume.total, 
		pow(PhysiCell::parameters.doubles("cell_radius"), 3.0) * 3.14159 * (4.0/3.0) * 2.0 
	) > UniformRandom() * 0.1);
}

bool Custom_cell::waiting_to_remove(Cell* cell, Phenotype& phenotype, double dt) {
	if (phenotype.cycle.data.elapsed_time_in_phase >= (8.6 * 60.0))
		return false;
	
	if (phenotype.volume.total < PhysiCell_constants::cell_removal_threshold_volume) 
		return false;

	return true;
}

void Custom_cell::set_initial_volume(Cell_Definition cell_def, float radius) {

	set_total_volume((4.0 / 3.0 * M_PI)*radius*radius*radius); 
	phenotype.volume.target_solid_nuclear = cell_def.phenotype.volume.target_solid_nuclear;
	phenotype.volume.target_solid_cytoplasmic = cell_def.phenotype.volume.target_solid_cytoplasmic;
	phenotype.volume.rupture_volume = cell_def.phenotype.volume.rupture_volume;
}

void Custom_cell::add_cell_basement_membrane_interactions( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);

	//Note that the distance_to_membrane function must set displacement values (as a normal vector)
	double max_interactive_distance = PhysiCell::parameters.doubles("max_interaction_factor") * phenotype.geometry.radius;
		
	double temp_a=0;
	double membrane_length = PhysiCell::parameters.doubles("membrane_length");
	double distance = pCustomCell->distance_to_membrane(pCell, phenotype, membrane_length);
	// Adhesion to basement membrane
	if(distance< max_interactive_distance)
	{
		temp_a= (1- distance/max_interactive_distance);
		temp_a*=temp_a;
		temp_a*=- PhysiCell::parameters.doubles("cell_basement_membrane_adhesion");
	}
	// Repulsion from basement membrane
	double temp_r=0;
	if ( distance < phenotype.geometry.radius )
	{
		temp_r= (1- distance/phenotype.geometry.radius);
		temp_r*=temp_r;
		temp_r*= PhysiCell::parameters.doubles("cell_basement_membrane_repulsion");
	}
	temp_r+=temp_a;
	if(temp_r==0)
		return;

	pCell->velocity += temp_r * pCell->displacement;	
	return;	
}

/* Calculate agent distance to BM if defined */
double Custom_cell::distance_to_membrane(Cell* pCell, Phenotype& phenotype, double l)
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);
	std::string shape = pCustomCell->get_shape();
	if ( l > 0 )
	{
		if ( shape == "duct" )
			return pCustomCell->distance_to_membrane_duct(l);
		else if ( shape == "sphere" )
			return pCustomCell->distance_to_membrane_sphere(l);
		else if ( shape == "sheet" )
			return pCustomCell->distance_to_membrane_sheet(l);
	}
	return 0;
}


/* Distance to membrane Sphere 
 * Basement membrane is a sphere of radius BM_radius 
 * Sphere center is (0,0,0)
 * */
double Custom_cell::distance_to_membrane_sphere(double length)
{
	double distance_to_origin = sqrt( position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);  // distance to the origin 
	distance_to_origin = std::max(distance_to_origin, EPSILON);	  // prevents division by zero
	displacement = -1 / distance_to_origin * position;
	if ( (length - distance_to_origin) < 0 )
		displacement *= 2.0; // penalize more outside of the sphere cells, stronger rappel
	return fabs(length - distance_to_origin);
}

/* Distance to membrane Sheet
 * Basement membrane is a sheet of height 2*BM_radius 
 * Z value is in between -BM_radius and +BM_radius
 * */
double Custom_cell::distance_to_membrane_sheet(double length)
{
	double distance = fabs(position[2]);  // |z| position
	distance = std::max(distance, EPSILON);	  // prevents division by zero
	displacement[0] = 0;
	displacement[1] = 0;
	displacement[2] = -1 / distance * position[2];
	if ( (length - distance) < 0 )
		displacement *= 2.0; // penalize more outside of the sphere cells, stronger rappel
	return fabs(length - distance);
}

/// Distance to membrane functions
double Custom_cell::distance_to_membrane_duct(double length)
{
	//Note that this function assumes that duct cap center is located at <0, 0, 0>
	if ( position[0] >= 0 ) // Cell is within the cylinder part of the duct
	{
		double distance_to_x_axis= sqrt( position[1]*position[1] + position[2]*position[2]);
		distance_to_x_axis = std::max(distance_to_x_axis, EPSILON);		// prevents division by zero
		displacement[0]=0; 
		displacement[1]= -position[1]/ distance_to_x_axis; 
		displacement[2]= -position[2]/ distance_to_x_axis; 
		return fabs(length - distance_to_x_axis);
	}

	// Cell is inside the cap of the duct
	double distance_to_origin= sqrt( position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);  // distance to the origin 
	distance_to_origin = std::max(distance_to_origin, EPSILON);			  // prevents division by zero
	displacement = -1 / distance_to_origin * position;
	return fabs(length - distance_to_origin);
}

/* Return if level of protein given by index around the cell is high enough (compared to given threshold) */
int Custom_cell::feel_enough(std::string field, Custom_cell pCell)
{	
	int voxel_index = get_current_mechanics_voxel_index();
	//return local_density(field) > cell_line->prot_threshold; 
	int ind = BioFVM::microenvironment.find_density_index( field );
	if ( ind >= 0 )	
		return BioFVM::microenvironment.density_vector(voxel_index)[ind] > get_threshold(field); 
	return -1;
}

/* Return true if level of oxygen is lower than necrosis critical level */
bool Custom_cell::necrotic_oxygen()
{	
	int oxygen_substrate_index = BioFVM::microenvironment.find_density_index( "oxygen" );
	double ox = this->nearest_density_vector()[oxygen_substrate_index];
	//std::cout << ox << " " << (this->parameters.o2_necrosis_threshold) - ox << std::endl;
	if ( ox >= 0 )	
		return ( UniformRandom() * 5 < (this->parameters.o2_necrosis_threshold - ox) );
   return false;	
}

void Custom_cell::add_TGFbeta_interaction(int index_voxel){
	int TGFbeta_index = BioFVM::microenvironment.find_density_index( "TGFbeta" );
	int ecm_index = BioFVM::microenvironment.find_density_index( "ecm" );
	double dens_tgfb = get_microenvironment()->nearest_density_vector(index_voxel)[TGFbeta_index];
	double dens_ecm = get_microenvironment()->nearest_density_vector(index_voxel)[ecm_index];
	double ratio = dens_ecm / dens_tgfb;
	if(ratio < PhysiCell::parameters.doubles("ECM_TGFbeta_ratio")){
	//std::cout << "dens =" << dens << std::endl;
	dens_tgfb = std::min( dens_tgfb, 1.0 );
	//std::cout << "dens =" << dens << std::endl; 
	if ( dens_tgfb > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		displacement = position - get_container()->underlying_mesh.voxels[index_voxel].center;
		double distance = norm(displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);

		double max_interactive_distance = (PhysiCell::parameters.doubles("max_interaction_factor")*phenotype.geometry.radius) + ecmrad;
		if ( distance < max_interactive_distance ) 
		{	
			TGFbeta_contact += dens_tgfb * (max_interactive_distance-distance);
			//std::cout << "TGFbeta =" << TGFbeta_contact << std::endl; 
		}
	}
	}
}