#include "custom_cell.h"

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


