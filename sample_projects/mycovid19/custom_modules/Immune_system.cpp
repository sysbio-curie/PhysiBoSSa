#include "./Immune_system.h" 
#include "Macrophage.h"
#include "TCell.h"

using namespace PhysiCell;

void Immune_system::create_infiltrating_Tcell(void)
{
	// static Cell_Definition* pCD = find_cell_definition( "CD8 Tcell" );
	// create_infiltrating_immune_cell( pCD ); 

	// return; 
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static bool setup_done = false; 
	
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// keep cells away from the outer edge 
	
	if( setup_done == false )
	{
		Xmin += 0.01*Xrange; 
		Ymin += 0.01*Yrange; 
		Zmin = 0;
		
		Xrange *= 0.98;
		Yrange *= 0.98;
		Zrange = 0.0; 
		setup_done = true; 
	}
	
	// create some of each type of cell 
	
	Cell* pC;
	
	std::vector<double> position = {0,0,0}; 
	position[0] = Xmin + UniformRandom()*Xrange; 
	position[1] = Ymin + UniformRandom()*Yrange; 
	position[2] = parameters.doubles("immune_z_offset"); 
		
	pC = create_cell( get_cell_definition("CD8 Tcell" ) ); 
	pC->assign_position( position );
	pC->update_voxel_in_container();

	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	pC->phenotype.secretion.uptake_rates[proinflammatory_cytokine_index] = 
		parameters.doubles("activated_cell_cytokine_uptake_rate"); // 10;
	
	return;
}

void Immune_system::immune_cell_recruitment( double dt )
{
	static int proinflammatory_cytokine_index = 
		microenvironment.find_density_index("pro-inflammatory cytokine");
	
	static double dt_immune = parameters.doubles( "immune_dt" ); 
	static double t_immune = 0.0; 
	static double t_last_immune = 0.0; 
	static double t_next_immune = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_immune > t_next_immune- tolerance )
	{
		double elapsed_time = (t_immune - t_last_immune );
//		std::cout << "Immune time! " << t_immune << " (elapsed: " << elapsed_time << ") " << std::endl; 
		
		
		// CD8 T cell recruitment 
		
		static double CD8_Tcell_recruitment_rate = parameters.doubles( "CD8_Tcell_max_recruitment_rate" ); 
		static double TC_min_signal = parameters.doubles( "CD8_Tcell_recruitment_min_signal" ); 
		static double TC_sat_signal = parameters.doubles( "CD8_Tcell_recruitment_saturation_signal" ); 
		static double TC_max_minus_min = TC_sat_signal - TC_min_signal; 
		
		double total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
		double total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[proinflammatory_cytokine_index] - TC_min_signal ); 
			dRate /= TC_max_minus_min; 
			// crop to [0,1] 
			if( dRate > 1 ) 
			{ dRate = 1; } 
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate; 
		}	
		// multiply by dV and rate_max 
		total_scaled_signal = total_rate; 
		
		total_rate *= microenvironment.mesh.dV; 
		total_rate *= CD8_Tcell_recruitment_rate; 
		
		// expected number of new neutrophils 
		int number_of_new_cells = (int) round( total_rate * elapsed_time ); 
		if( number_of_new_cells )
		{
			// std::cout << "\tRecruiting " << number_of_new_cells << " CD8 T cells ... " << std::endl; 
			
//			std::cout << "\tTotal signal/dV : " << total_scaled_signal << std::endl;
//			std::cout << "\tTotal signa : " << total_scaled_signal * microenvironment.mesh.dV << std::endl; 
//			double total_volume = microenvironment.mesh.dV * microenvironment.mesh.voxels.size() ; 
//			std::cout << "\tmean signal : " << total_scaled_signal * microenvironment.mesh.dV / total_volume << std::endl; 
			for( int n = 0; n < number_of_new_cells ; n++ )
			{ create_infiltrating_Tcell(); }
		}
		
		t_last_immune = t_immune; 
		t_next_immune = t_immune + dt_immune; 
		
		// std::cout << "\t\tnext immune time: " << t_next_immune << std::endl;  
	}
	t_immune += dt; 
	return; 
}

void Immune_system::initialize() {
	
	Cell_Definition* macrophagesCD = find_cell_definition("macrophage"); 
	Macrophage::setup_cell_definition(macrophagesCD);
	macrophagesCD->phenotype.sync_to_microenvironment( &microenvironment );
	
	Cell_Definition* tcellCD = find_cell_definition("CD8 Tcell"); 
	TCell::setup_cell_definition(tcellCD);
	tcellCD->phenotype.sync_to_microenvironment( &microenvironment );
}

void Immune_system::create(){
	
	Cell* pC;

	static double x_min = microenvironment.mesh.bounding_box[0]; 
	static double y_min = microenvironment.mesh.bounding_box[1]; 

	static double x_max = microenvironment.mesh.bounding_box[3]; 
	static double y_max = microenvironment.mesh.bounding_box[4]; 

	// now place immune cells 
	for (int i=0; i < parameters.ints("number_macrophages"); i++) {
	
		pC = create_cell( get_cell_definition("macrophage" ) ); 
		pC->assign_position( 
			x_min + (x_max-x_min)*UniformRandom(),
			y_min + (y_max-y_min)*UniformRandom(),
			0.0 + parameters.doubles("immune_z_offset")
		);
		
	}	
}
	
	

void Immune_system::keep_immune_cells_off_edge( void )
{
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 

	static bool setup_done = false; 
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// warning hardcoded
	static double relative_edge_margin = 0; // 0.1; 
	static double relative_interior = 1 - 2 * relative_edge_margin; 
	
	if( setup_done == false )
	{
		Xmin += relative_edge_margin*Xrange; 
		Ymin += relative_edge_margin*Yrange; 
		Zmin += relative_edge_margin*Zrange;
		
		Xrange *= relative_interior;
		Yrange *= relative_interior;
		Zrange *= relative_interior;  
		setup_done = true; 
	}
	
	static int epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	
	for( int n=0 ; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if (pC->type != epithelial_type)
			static_cast<Immune_cell*>(pC)->keep_in_bounds(Xmin, Xrange, Ymin, Yrange, Zmin, Zrange);
		else{
			if (pC->position[2] != 0)	std::cout << "We are not on zero :(" << std::endl;
		}
	}
	return; 
}


void Immune_system::keep_immune_cells_in_bounds( double dt )
{
	static double dt_bounds = 0.1; 
	static double next_time = 0.0; 

	static double t_bounds = 0.0; 
	static double t_last_bounds = 0.0; 
	static double t_next_bounds = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_bounds > t_next_bounds- tolerance )
	{
		double elapsed_time = (t_bounds - t_last_bounds );
		
		keep_immune_cells_off_edge(); 
		
		t_last_bounds = t_bounds; 
		t_next_bounds = t_bounds + dt_bounds; 
	}
	t_bounds += dt; 

	return; 
}
