#include "custom_cell.h"

Custom_cell::Custom_cell() {
    pintegrin = 0.5;
	pmotility = 0.5;
	padhesion = 0.5;
	mmped = 0;
	motility.resize(3, 0.0);
}

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void Custom_cell::add_ecm_interaction( int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	//if (dens > PhysiCell::EPSILON || dens2 > PhysiCell::EPSILON) { std::cout << dens << "    " << dens2 << std::endl;};
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > PhysiCell::EPSILON )
	{
		// Distance between agent center and ECM voxel center
		displacement = position - get_container()->underlying_mesh.voxels[index_voxel].center;
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
	Cecm[0] = PhysiCell::parameters.ints("ecm_adhesion_min");
	Cecm[1] = PhysiCell::parameters.ints("ecm_adhesion_min");
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
	Ccca_heterotypic[0] = PhysiCell::parameters.doubles("heterotypic_adhesion_min");
	Ccca_heterotypic[1] = PhysiCell::parameters.doubles("heterotypic_adhesion_max");
	Ccca_homotypic[0] = PhysiCell::parameters.doubles("homotypic_adhesion_min");
	Ccca_homotypic[1] = PhysiCell::parameters.doubles("homotypic_adhesion_max");
	double adh = 0;
	if ( &(Cell::phenotype) == &(other_cell->Cell::phenotype) )
		adh = std::min( get_homotypic_strength(padhesion), custom_other_cell->get_homotypic_strength(padhesion) );
	else
		adh = std::min( get_heterotypic_strength(padhesion), custom_other_cell->get_heterotypic_strength(padhesion) );

	return adh;
}


double Custom_cell::get_adhesion()
{
	return 1;
}


/* Motility with random direction, and magnitude of motion given by customed coefficient */
void Custom_cell::set_3D_random_motility( double dt )
{
    double probability = UniformRandom();
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
