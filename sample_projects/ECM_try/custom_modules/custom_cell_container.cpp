#include "../BioFVM/BioFVM_agent_container.h"
#include "../core/PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "custom_cell.h"
#include "custom_cell_container.h"
using namespace BioFVM;
using namespace PhysiCell;

std::vector<Custom_cell*> *all_custom_cells;

Custom_cell_container::Custom_cell_container()
{

all_custom_cells = (std::vector<Custom_cell*> *) &all_cells;
custom_cells_ready_to_divide = (std::vector<Custom_cell*> *) &cells_ready_to_divide;
custom_cells_ready_to_die= (std::vector<Custom_cell*> *) &cells_ready_to_die;

return;

};


Custom_cell_container* create_custom_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size )
{
	Custom_cell_container* custom_cell_container = new Custom_cell_container;
	custom_cell_container->initialize( m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
		m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
		m.mesh.bounding_box[2], m.mesh.bounding_box[5],  mechanics_voxel_size );
	m.agent_container = (Agent_Container*) custom_cell_container; 
	
	if( BioFVM::get_default_microenvironment() == NULL )
	{ 
		BioFVM::set_default_microenvironment( &m ); 
	}
	
	return custom_cell_container; 
}

void Custom_cell_container::update_all_cells(double t)
{
	// update_all_cells(t, dt_settings.cell_cycle_dt_default, dt_settings.mechanics_dt_default);
	
	update_all_cells(t, phenotype_dt, mechanics_dt , diffusion_dt );
	
	return; 
}

void Custom_cell_container::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ )
{
	// secretions and uptakes. Syncing with BioFVM is automated. 

	#pragma omp parallel for 
	for( int i=0; i < (*all_custom_cells).size(); i++ )
	{
		(*all_custom_cells)[i]->phenotype.secretion.advance( (*all_custom_cells)[i], (*all_custom_cells)[i]->phenotype , diffusion_dt_ );
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
		for( int i=0; i < (*all_custom_cells).size(); i++ )
		{
			if( (*all_custom_cells)[i]->is_out_of_domain == false )
			{
				(*all_custom_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle ); 
			}
		}
		
		// process divides / removes 
		for( int i=0; i < (*custom_cells_ready_to_divide).size(); i++ )
		{
			(*custom_cells_ready_to_divide)[i]->divide();
		}
		for( int i=0; i < (*custom_cells_ready_to_die).size(); i++ )
		{	
			(*custom_cells_ready_to_die)[i]->die();	
		}
		num_divisions_in_current_step+=  (*custom_cells_ready_to_divide).size();
		num_deaths_in_current_step+=  (*custom_cells_ready_to_die).size();
		
		(*custom_cells_ready_to_die).clear();
		(*custom_cells_ready_to_divide).clear();
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
		for( int i=0; i < (*all_custom_cells).size(); i++ )
		{

			if(!(*all_custom_cells)[i]->is_out_of_domain && (*all_custom_cells)[i]->is_movable && (*all_custom_cells)[i]->functions.update_velocity )
			{
				// update_velocity already includes the motility update 
				//(*all_cells)[i]->phenotype.motility.update_motility_vector( (*all_cells)[i] ,(*all_cells)[i]->phenotype , time_since_last_mechanics ); 
				(*all_custom_cells)[i]->functions.update_velocity( (*all_custom_cells)[i], (*all_custom_cells)[i]->phenotype, time_since_last_mechanics);
			}

			if ( !(*all_custom_cells)[i]->passive() )
				((*all_custom_cells)[i])->degrade_ecm( mechanics_dt_ );

			if( (*all_custom_cells)[i]->functions.custom_cell_rule )
			{
				(*all_custom_cells)[i]->functions.custom_cell_rule((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}
		}
		// Calculate new positions
		#pragma omp parallel for 
		for( int i=0; i < (*all_custom_cells).size(); i++ )
		{
			if(!(*all_custom_cells)[i]->is_out_of_domain && (*all_custom_cells)[i]->is_movable)
			{
				(*all_custom_cells)[i]->update_position(time_since_last_mechanics);
			}
		}
		
		// When somebody reviews this code, let's add proper braces for clarity!!! 
		
		// Update cell indices in the container
		for( int i=0; i < (*all_custom_cells).size(); i++ )
			if(!(*all_custom_cells)[i]->is_out_of_domain && (*all_custom_cells)[i]->is_movable)
				(*all_custom_cells)[i]->update_voxel_in_container();
		last_mechanics_time=t;
	}
	
	initialzed=true;
	return;
}
