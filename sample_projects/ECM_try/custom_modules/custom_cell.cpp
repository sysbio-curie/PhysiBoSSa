#include "custom_cell.h"
#include "custom_cell_container.h"
#include "../core/PhysiCell_utilities.h"
#include "../core/PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h" 
#include<limits.h>


using namespace BioFVM; 


Custom_cell::Custom_cell()
{
	displacement.resize(3, 0.0); // state? 
	pintegrin = 0.5;
	pmotility = 0.5;
	padhesion = 0.5;
	mmped = 0;
	return; 
}

Custom_cell_container * Custom_cell::get_custom_container()
{
	if(custom_container == NULL)
	{
		custom_container = (Custom_cell_container *)get_microenvironment()->agent_container;
	}
	
	return custom_container;
}

void Custom_cell::advance_bundled_phenotype_functions( double dt_ )
{
	// call the custom code to update the phenotype 
	if( functions.update_phenotype )
	{	functions.update_phenotype( this , phenotype , dt_ ); }
	
	// update volume 
	if( functions.volume_update_function )
	{
		functions.volume_update_function(this,phenotype,dt_); 
		
		// The following line is needed in every volume 
		// regulation method (it sets BioFVM total_volume)
		
		set_total_volume( phenotype.volume.total ); 
	}
	
	// update geometry
	phenotype.geometry.update( this, phenotype, dt_ );
	
	// check for new death events 
	if( phenotype.death.check_for_death( dt_ ) == true )
	{
		// if so, change the cycle model to the current death model 
		phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() ); 
		
		// also, turn off motility.
		
		phenotype.motility.is_motile = false; 
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		functions.update_migration_bias = NULL;
		
		// turn off secretion, and reduce uptake by a factor of 10 
		phenotype.secretion.set_all_secretion_to_zero();
		phenotype.secretion.scale_all_uptake_by_factor( 0.10 );
		
		// make sure to run the death entry function 
		if( phenotype.cycle.current_phase().entry_function )
		{
			phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ ); 
		}
	}
	
	// advance cycle model (for both cell cycle and death cycle models)
	phenotype.cycle.advance_cycle( this, phenotype, dt_ ); 
	if( phenotype.flagged_for_removal )
	{
		custom_flag_for_removal(); 
		phenotype.flagged_for_removal = false; 
	}
	if( phenotype.flagged_for_division )
	{
		custom_flag_for_division(); 
		phenotype.flagged_for_division = false; 
	}
	
	return; 
}


void Custom_cell::custom_flag_for_division( void )
{
	get_custom_container()->custom_flag_cell_for_division( this );
	return; 
}

void Custom_cell::custom_flag_for_removal( void )
{
	get_custom_container()->custom_flag_cell_for_removal( this );
	return;
}


void Custom_cell::add_potentials(Custom_cell* other_agent)
{
	if( this->ID == other_agent->ID )
	{ return; }

	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		displacement[i] = position[i] - (*other_agent).position[i]; 
		distance += displacement[i] * displacement[i]; 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = phenotype.geometry.radius+ (*other_agent).phenotype.geometry.radius; 
	
	double RN = phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	// else if( distance < RN ) 
	// {
		// double M = 1.0; 
		// c = 1.0 - RN/R; 
		// c *= c; 
		// c -= M; 
		// temp_r = ( c*distance/RN  + M  ); 
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
	temp_r *= effective_repulsion; 
	
	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original 
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius + 
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;
		
	double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance * phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
		// hyp: junction strength is limited by weakest cell
		double adh;
		double thisadh = get_adhesion();
		double otadh = other_agent->get_adhesion();

		// first case, passive cell with active cell
		if ( thisadh == 0 && otadh == 1 )
		{
			ecm_contact += (max_interactive_distance-distance);
			adh = static_cast<Cell*>(other_agent)->integrinStrength();
		}
		else
		{
			// second case, active cell with passive cell
			if ( thisadh == 1 && otadh == 0 )
			{
				ecm_contact += (max_interactive_distance-distance);
				adh = static_cast<Cell*>(this)->integrinStrength();
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
					cell_contact += (max_interactive_distance-distance);
					adh = ( static_cast<Cell*>(this) )->adhesion( static_cast<Cell*>(other_agent) );
				}
			}
		}
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ ) 
	// {
	//	velocity[i] += displacement[i] * temp_r; 
	// }
	axpy( &velocity , temp_r , displacement ); 
	
	return;
}

Custom_cell* create_custom_cell( void )
{
	Custom_cell* pNew; 
	pNew = new Custom_cell;		
	(*all_cells).push_back( pNew ); 
	pNew->index=(*all_cells).size()-1;
	
	// new usability enhancements in May 2017 
	
	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}

	// All the phenotype and other data structures are already set 
	// by virtue of the default Cell constructor. 
	
	return pNew; 
}

// in the future, I might swap this with create_cell(): 
// In that "create_cell()" uses "create_cell( cell_defaults )" 
Custom_cell* create_custom_cell( Cell_Definition& cd )
{
	Custom_cell* pNew = create_custom_cell(); 
	
	// use the cell defaults; 
	pNew->type = cd.type; 
	pNew->type_name = cd.name; 
	
	pNew->custom_data = cd.custom_data; 
	pNew->parameters = cd.parameters; 
	pNew->functions = cd.functions; 
	
	pNew->phenotype = cd.phenotype; 
	pNew->is_movable = true;
	pNew->is_out_of_domain = false;
	pNew->displacement.resize(3,0.0); // state? 
	
	pNew->assign_orientation();
	
	return pNew; 
}



/* Degrade the surrounding ECM 
 *
 * param dt time step */
void Custom_cell::degrade_ecm( double dt )
{
	if ( is_out_of_domain )
		return;
	if ( !mmped ) 
		return;

	// Check if there is ECM material in given voxel
	int ecm_index = get_microenvironment()->find_density_index("ecm");
	int current_index = get_current_mechanics_voxel_index();
	#pragma omp critical
	{
		double dens = get_microenvironment()->nearest_density_vector(current_index)[ecm_index];
		if ( dens > EPSILON )
		{
			dens -= (PhysiCell::parameters.ints("ecm_degradation") * pintegrin) * dt; // to change by a rate
			dens = dens > 0 ? dens : 0;
			get_microenvironment()->nearest_density_vector(current_index)[ecm_index] = dens;
		}
	}
}


/* Return value of adhesion strength with ECM according to integrin level */
double Custom_cell::integrinStrength()
{ 
	Cecm[0] = PhysiCell::parameters.ints("ecm_adhesion_min");
	Cecm[1] = PhysiCell::parameters.ints("ecm_adhesion_min");
	return Cell::get_integrin_strength( pintegrin ); 
}


/* Return if cell has enough contact with other cells (compared to given threshold determined by the given level) */	
bool Custom_cell::has_neighbor(int level)
{ 
	if ( level == 0 )
		return contact_cell() > PhysiCell::parameters.doubles("contact_cell_cell_threshold"); 
	else
		return contact_cell() > (2 * PhysiCell::parameters.doubles("contact_cell_cell_threshold")); 
}


/* Return level of protein given by index around the cell */
double Custom_cell::local_density(std::string field)
{ 
	int ind = BioFVM::microenvironment.find_density_index(field);
	if ( ind >= 0 )
		return (nearest_density_vector())[ind]; 
	return -1;
}


/* Return true if level of oxygen is lower than necrosis critical level */
bool Custom_cell::necrotic_oxygen()
{
	double ox = local_density("oxygen");
	//std::cout << ox << " " << (cell_line->o2_necrotic) - ox << std::endl;
	if ( ox >= 0 )	
		return ( PhysiCell::UniformRandom() * 0.005 < (PhysiCell::parameters.ints("oxygen_necrotic") - ox) );
   return false;	
}

/* Motility with random direction, and magnitude of motion given by customed coefficient */
void Custom_cell::set_3D_random_motility( double dt )
{
	double probability = PhysiCell::UniformRandom();
	motility_magnitude[0] = PhysiCell::parameters.doubles("motility_amplitude_min");
	motility_magnitude[1] = PhysiCell::parameters.doubles("motility_amplitude_max");
	if ( probability < dt / PhysiCell::parameters.doubles("persistence") )
	{
		std::vector<double> tmp;
		double temp_angle = 2 * M_PI * PhysiCell::UniformRandom();
		double temp_phi = M_PI * PhysiCell::UniformRandom();
		tmp[0] = cos( temp_angle ) * sin( temp_phi );
		tmp[1] = sin( temp_angle ) * sin( temp_phi );
		tmp[2] = cos( temp_phi );
		motility = get_motility_amplitude(pmotility) * tmp;
	}
}

/*
 * Motility in the polarity axis migration + little noise
 * Persistence in the update polarization
 * */
void Custom_cell::set_3D_polarized_motility( double dt )
{
	// mot = (1-p) * r + p * pol
	double temp_angle = 2 * M_PI * PhysiCell::UniformRandom();
	double temp_phi = M_PI * PhysiCell::UniformRandom();
	motility[0] = cos( temp_angle ) * sin( temp_phi );
	motility[1] = sin( temp_angle ) * sin( temp_phi );
	motility[2] = cos( temp_phi );
	motility *= (1 - PhysiCell::parameters.doubles("polarity_coefficient"));
	std::vector<double> polarization;
	polarization[0]= state.orientation[0];
	polarization[1]= state.orientation[1];
	polarization[2]= state.orientation[2];
	std::vector<double> pol_dir;
	double pol_norm = norm(polarization); //normal to polaization used to calculate the vestor direction for polarization
	pol_dir[0] = polarization[0]/pol_norm;
	pol_dir[1] = polarization[1]/pol_norm;
	pol_dir[2] = polarization[2]/pol_norm;
	motility += PhysiCell::parameters.doubles("polarity_coefficient") * pol_dir;
	// Normalized it
	normalize(motility);
	// mot = mot_coef * mot_dir
	motility *= get_motility_amplitude(pmotility);
}



/**
 * Calculate motility forces according to mode:
 * 0, random; 1, along polarity axis; other: nothing
 * */
void Custom_cell::set_motility( double dt )
{
	// Cell frozen, cannot actively move
	if ( freezed > 2 )
		return;
	switch( PhysiCell::parameters.ints("mode_motility") )
	{
		case 0:
			set_3D_random_motility(dt);
			break;
		case 1:
			set_3D_polarized_motility(dt);
			break;
		default:
			return;
			break;
	}
	velocity += motility;
}

/* Update the value of freezing of the cell with bitwise operation
* Do a bitwise-or comparison on freezed and input parameter:
* if freezed = 0, it will be the value of the parameter frozen
* if freezed = 1, it will be either 1 (frozen = 0) or 3 (frozen = 3) */
void Custom_cell::freezer( int frozen )
{
	freezed = freezed | frozen;
}

double Custom_cell::get_adhesion()
{
	return 1;
}

/* Calculate adhesion coefficient with other cell */
double Custom_cell::adhesion( Cell* other_cell )
{
	Ccca_heterotypic[0] = PhysiCell::parameters.doubles("heterotypic_adhesion_min");
	Ccca_heterotypic[1] = PhysiCell::parameters.doubles("heterotypic_adhesion_max");
	Ccca_homotypic[0] = PhysiCell::parameters.doubles("homotypic_adhesion_min");
	Ccca_homotypic[1] = PhysiCell::parameters.doubles("homotypic_adhesion_max");
	double adh = 0;
	if ( &(Cell::phenotype) == &(other_cell->Cell::phenotype) )
		adh = std::min( get_homotypic_strength(padhesion), other_cell->get_homotypic_strength(padhesion) );
	else
		adh = std::min( get_heterotypic_strength(padhesion), other_cell->get_heterotypic_strength(padhesion) );

	return adh;
}


void Custom_cell::add_cell_basement_membrane_interactions( double dt, double distance ) 
{
	double rad = phenotype.geometry.radius;
	//Note that the distance_to_membrane function must set displacement values (as a normal vector)
	double max_interactive_distance = PhysiCell::parameters.doubles("max_interaction_factor") * rad;
		
	double temp_a=0;
	// Adhesion to basement membrane
	if(distance< max_interactive_distance)
	{
		temp_a= (1- distance/max_interactive_distance);
		temp_a*=temp_a;
		temp_a*=- PhysiCell::parameters.doubles("cell_basement_membrane_adhesion");
	}
	// Repulsion from basement membrane
	double temp_r=0;
	if ( distance < rad )
	{
		temp_r= (1- distance/rad);
		temp_r*=temp_r;
		temp_r*= PhysiCell::parameters.doubles("cell_basement_membrane_repulsion");
	}
	temp_r+=temp_a;
	if(temp_r==0)
		return;

	velocity += temp_r * displacement;	
	return;	
}


/// Distance to membrane functions
double Custom_cell::distance_to_membrane_duct(double length)
{
	//Note that this function assumes that duct cap center is located at <0, 0, 0>
	if ( position[0] >= 0 ) // Cell is within the cylinder part of the duct
	{
		double distance_to_x_axis= sqrt((position[1] * position[1]) + (position[2] * position[2]));
		distance_to_x_axis = std::max(distance_to_x_axis, EPSILON);		// prevents division by zero
		displacement[0]=0; 
		displacement[1]= -position[1]/ distance_to_x_axis; 
		displacement[2]= -position[2]/ distance_to_x_axis; 
		return fabs(length - distance_to_x_axis);
	}

	// Cell is inside the cap of the duct
	double distance_to_origin= norm(position);  // distance to the origin 
	distance_to_origin = std::max(distance_to_origin, EPSILON);			  // prevents division by zero
	displacement = -1 / distance_to_origin * position;
	return fabs(length - distance_to_origin);
}

/* Distance to membrane Sphere 
 * Basement membrane is a sphere of radius BM_radius 
 * Sphere center is (0,0,0)
 * */
double Custom_cell::distance_to_membrane_sphere(double length)
{
	double distance_to_origin = norm(position);  // distance to the origin 
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

/* Calculate agent distance to BM if defined */
double Custom_cell::distance_to_membrane(double l, std::string shape)
{
	if ( l > 0 )
	{
		if ( shape == "duct" )
			return distance_to_membrane_duct(l);
		else if ( shape == "sphere" )
			return distance_to_membrane_sphere(l);
		else if ( shape == "sheet" )
			return distance_to_membrane_sheet(l);
	}
	return 0;
}

/* Update cell velocity */
void Custom_cell::update_cell_motion( double time_since_last, double l, std::string shape )
{
	cell_contact = 0;
	ecm_contact = 0;
	nucleus_deform = 0;
	if( !is_out_of_domain && is_movable)
		update_velocity( time_since_last, l, shape );
}


void Custom_cell::update_velocity( double dt, double l, std::string shape ) 
{
	double dist = distance_to_membrane(l, shape);
	if ( dist > 0 )
		add_cell_basement_membrane_interactions(dt, dist);


	//First check the neighbors in my current voxel
	for( Custom_cell * neighbor : get_custom_container()->agent_grid[get_current_mechanics_voxel_index()] )
	{
		add_potentials( neighbor );
	}
	
	int ecm_index = BioFVM::microenvironment.find_density_index("ecm");
	if ( ecm_index >= 0 )
		add_ecm_interaction( ecm_index, get_current_mechanics_voxel_index() );

	for ( auto neighbor_voxel_index : get_container()-> underlying_mesh.moore_connected_voxel_indices[get_current_mechanics_voxel_index()] )
	{
		if ( !PhysiCell::is_neighbor_voxel(this, (get_container()->underlying_mesh.voxels[get_current_mechanics_voxel_index()].center), (get_container()->underlying_mesh.voxels[neighbor_voxel_index].center), neighbor_voxel_index) )
			continue;
		if ( ecm_index >= 0 )
			add_ecm_interaction( ecm_index, neighbor_voxel_index );
		for( Custom_cell * other_neighbor : get_custom_container()->agent_grid[neighbor_voxel_index] )
		{
			add_potentials( other_neighbor );
		}
	}

	// Add active motility term
	if ( !passive() )
		set_motility(dt);
}


/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void Custom_cell::add_ecm_interaction( int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	double dens = get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > PhysiCell::EPSILON )
	{
		// Distance between agent center and ECM voxel center
		displacement = position - get_microenvironment()->get_voxel_center(index_voxel);
		double distance = norm(displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, PhysiCell::EPSILON);
		
		double ecmrad = get_microenvironment()->voxel_rad(index_voxel);
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
			temp_a *= dens * ( static_cast<Cell*>(this) )->integrinStrength();
			tmp_r -= temp_a;
		}
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		velocity += tmp_r * displacement;
	}
}

