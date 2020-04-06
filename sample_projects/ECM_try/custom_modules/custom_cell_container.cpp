#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"
#include "custom_cell_container.h"
using namespace BioFVM;
using namespace PhysiCell;

void Custom_cell_container::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ )
{
	// secretions and uptakes. Syncing with BioFVM is automated. 

	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		(*all_cells)[i]->phenotype.secretion.advance( (*all_cells)[i], (*all_cells)[i]->phenotype , diffusion_dt_ );
	}
	
	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	static double phenotype_dt_tolerance = 0.001 * phenotype_dt_; 
	static double mechanics_dt_tolerance = 0.001 * mechanics_dt_; 
	
	if( fabs(time_since_last_cycle-phenotype_dt_ ) < phenotype_dt_tolerance || !initialzed)
	{
		// Reset the max_radius in each voxel. It will be filled in set_total_volume
		// It might be better if we calculate it before mechanics each time 
		std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);
		
		if(!initialzed)
		{
			time_since_last_cycle = phenotype_dt_;
		}
		
		// new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update, 
		// checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if( (*all_cells)[i]->is_out_of_domain == false )
			{
				(*all_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle ); 
			}
		}
		
		// process divides / removes 
		for( int i=0; i < cells_ready_to_divide.size(); i++ )
		{
			cells_ready_to_divide[i]->divide();
		}
		for( int i=0; i < cells_ready_to_die.size(); i++ )
		{	
			cells_ready_to_die[i]->die();	
		}
		num_divisions_in_current_step+=  cells_ready_to_divide.size();
		num_deaths_in_current_step+=  cells_ready_to_die.size();
		
		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
		last_cell_cycle_time= t;
	}
		
	double time_since_last_mechanics= t- last_mechanics_time;
	
	// if( time_since_last_mechanics>= mechanics_dt || !initialzed)
	if( fabs(time_since_last_mechanics - mechanics_dt_) < mechanics_dt_tolerance || !initialzed)
	{
		if(!initialzed)
		{
			time_since_last_mechanics = mechanics_dt_;
		}
		
		// new February 2018 
		// if we need gradients, compute them
		if( default_microenvironment_options.calculate_gradients ) 
		{ microenvironment.compute_all_gradient_vectors();  }
		// end of new in Feb 2018 		
		
		// Compute velocities
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{

			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable && (*all_cells)[i]->functions.update_velocity )
			{
				// update_velocity already includes the motility update 
				//(*all_cells)[i]->phenotype.motility.update_motility_vector( (*all_cells)[i] ,(*all_cells)[i]->phenotype , time_since_last_mechanics ); 
				(*all_cells)[i]->functions.update_velocity( (*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}

			if ( !(*all_cells)[i]->passive() )
				((*all_cells)[i])->degrade_ecm( mechanics_dt_ );

			if( (*all_cells)[i]->functions.custom_cell_rule )
			{
				(*all_cells)[i]->functions.custom_cell_rule((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}
		}
		// Calculate new positions
		#pragma omp parallel for 
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
			{
				(*all_cells)[i]->update_position(time_since_last_mechanics);
			}
		}
		
		// When somebody reviews this code, let's add proper braces for clarity!!! 
		
		// Update cell indices in the container
		for( int i=0; i < (*all_cells).size(); i++ )
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
				(*all_cells)[i]->update_voxel_in_container();
		last_mechanics_time=t;
	}
	
	initialzed=true;
	return;
}

void Custom_cell_container::custom_flag_cell_for_division( Custom_cell* pCell )
{ 
	#pragma omp critical 
	{cells_ready_to_divide.push_back( pCell );} 
	return; 
}

void Custom_cell_container::custom_flag_cell_for_removal( Custom_cell* pCell )
{ 
	#pragma omp critical 
	{cells_ready_to_die.push_back( pCell );} 
	return; 
}