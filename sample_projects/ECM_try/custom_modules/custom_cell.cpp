#include "custom_cell.h"

Custom_cell::Custom_cell() {
	ecm_contact = 0;
	nucleus_deform = 0;
	freezed = 0;
    pintegrin = 0.5;
	pmotility = 0.5;
	padhesion = 0.5;
	mmped = 0;
	ecmrad = sqrt(3.0) * get_microenvironment()->mesh.dx / 2.0;
	motility.resize(3, 0.0);	
}

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void Custom_cell::add_ecm_interaction( int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
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
	return get_integrin_strength( pintegrin ); 
}

/* Return if cell has enough contact with other cells (compared to given threshold determined by the given level) */	
bool Custom_cell::has_neighbor(int level)
{ 
	if ( level == 0 )
		return contact_cell() > PhysiCell::parameters.doubles("contact_cell_cell_threshold"); 
	else
		return contact_cell() > (2 * PhysiCell::parameters.doubles("contact_cell_cell_threshold")); 
}

/* Calculate adhesion coefficient with other cell */
double Custom_cell::adhesion( Cell* other_cell )
{
    Custom_cell* custom_other_cell = static_cast<Custom_cell*>(other_cell);
	
	double adh = 0;
	if ( this->type == other_cell->type )
		adh = std::min( get_homotypic_strength(padhesion), custom_other_cell->get_homotypic_strength(padhesion) );
	else
		adh = std::min( get_heterotypic_strength(padhesion), custom_other_cell->get_heterotypic_strength(padhesion) );
	
	return adh;
}


double Custom_cell::get_adhesion()
{
	if (passive())
		return 0;
	else
		return 1;
}


/* Motility with random direction, and magnitude of motion given by customed coefficient */
void Custom_cell::set_3D_random_motility( double dt )
{
    double probability = UniformRandom();
	
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
    polarization.resize(3, 0.0);
    polarization[0]= state.orientation[0];
    polarization[1]= state.orientation[1];
    polarization[2]= state.orientation[2];
    std::vector<double> pol_dir;
    pol_dir.resize(3, 0.0);
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

Cell* Custom_cell::create_custom_cell()
{
	Custom_cell* pNew; 
	pNew = new Custom_cell;		
	return static_cast<Cell*>(pNew); 
}

// Here I'm hoping that the argument used, time_since_last_mechanics, has the same value
// as mechanics_dt_. I should probably check later...
void Custom_cell::check_passive(Cell* cell, Phenotype& phenotype, double dt) {
	Custom_cell* t_cell = static_cast<Custom_cell*>(cell);
	if (!(t_cell->passive())) {
		t_cell->degrade_ecm(dt);
	}
}

void Custom_cell::custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	Custom_cell* pCustomCell = static_cast<Custom_cell*>(pCell);

	pCustomCell->ecm_contact = 0;
	pCustomCell->nucleus_deform = 0;

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

	for (auto neighbor_voxel_index : pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()])
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[neighbor_voxel_index].center, neighbor_voxel_index))
			continue;

		if ( ecm_index >= 0 )
			pCustomCell->add_ecm_interaction( ecm_index, neighbor_voxel_index );
	
		for( auto other_neighbor : pCell->get_container()->agent_grid[neighbor_voxel_index] )
		{
			pCell->add_potentials(other_neighbor);
		}
	}
	
	// Add active motility term
	if ( !(pCustomCell->passive()) )
		pCustomCell->set_motility(dt);
	
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
				custom_pCell->cell_contact += distance;
				adh = custom_pCell->adhesion(custom_otherCell);
			}
		}
	}
	return adh;
}

bool Custom_cell::wait_for_nucleus_growth (Cell* cell, Phenotype& phenotype, double dt) {
	return relative_diff( 
		phenotype.volume.total, 
		pow(PhysiCell::parameters.doubles("cell_radius"), 3.0) * 3.14159 * (4.0/3.0) * 2.0 
	) > UniformRandom() * 0.1;
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