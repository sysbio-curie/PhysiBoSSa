

#include "./tnf_response_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information tnf_dynamics_info;

void tnf_dynamics_model_setup()
{
    tnf_dynamics_info.name = "Tumor Necrotic Factor model dynamics"; 
	tnf_dynamics_info.version = "0.1.0";
	
    tnf_dynamics_info.main_function= tnf_dynamics_model; 

	// what custom data do I need?
	tnf_dynamics_info.cell_variables.push_back( "TNFR_activation_threshold" );

    tnf_dynamics_info.cell_variables.push_back( "unbound_external_TNFR" );
	tnf_dynamics_info.cell_variables.push_back( "bound_external_TNFR" );
	tnf_dynamics_info.cell_variables.push_back( "bound_internal_TNFR" );
    
    tnf_dynamics_info.cell_variables.push_back( "TNFR_binding_rate" ); 
	tnf_dynamics_info.cell_variables.push_back( "TNFR_endocytosis_rate" );
    tnf_dynamics_info.cell_variables.push_back( "TNFR_recycling_rate" );
	tnf_dynamics_info.cell_variables.push_back( "TFN_net_production_rate" );

	tnf_dynamics_info.register_model();
}


void tnf_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nTNF_external = microenvironment.find_density_index( "tnf" );

    static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_TNFR" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_TNFR" );
	static int nR_IB = pCell->custom_data.find_variable_index( "bound_internal_TNFR" );

    static int nR_bind = pCell->custom_data.find_variable_index( "TNFR_binding_rate" );
	static double R_binding_rate = pCell->custom_data[nR_bind];

	static int nR_endo = pCell->custom_data.find_variable_index( "TNFR_endocytosis_rate" ); 
	static double R_endo_rate = pCell->custom_data[nR_endo];

	static int nR_recycle = pCell->custom_data.find_variable_index( "TNFR_recycling_rate" ); 
	static double R_recyc_rate = pCell->custom_data[nR_recycle];
	

	if( phenotype.death.dead == true )
	{ return; } 
		
    // internalized TNF tells us how many have recently bound to receptors
	// TNF is internalized at:
	// phenotype.secretion.uptake_rates[nTNF_external] = 
	// 					pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU];
	
	// The internalization is only used to track the TNF
	// The following part of the code takes care of correcly managed
	double dR_EB = phenotype.molecular.internalized_total_substrates[nTNF_external];
	
	// if it tries to bind more TNF than there are receptors, compensate 
	if( dR_EB > pCell->custom_data[nR_EU] )
	{
		double excess_binding = dR_EB - pCell->custom_data[nR_EU];
		dR_EB = pCell->custom_data[nR_EU]; 
		// dump any excess back into the microenvironment
		static double to_density = 1.0 / microenvironment.mesh.dV; 
		// this needs omp critical because 2 cells writing to 1 voxel is not thread safe 
		#pragma omp critical 
		{ pCell->nearest_density_vector()[nTNF_external] += excess_binding * to_density; }
	}

	// Remove all the internalized TNF from cell
	phenotype.molecular.internalized_total_substrates[nTNF_external] = 0.0; 
	
	// Endocytosis 
	// The bounded receptor is internalized at a rate R_endo_rate
	double dR_IB = dt * R_endo_rate  * pCell->custom_data[nR_EB];
    if( dR_IB > pCell->custom_data[nR_EB] )
	{ dR_IB = pCell->custom_data[nR_EB]; }

	// Recylcing
	// The internalized bounded TNFR release the TNF
    // The TNF is instantaneously degraded by the cell
    // The TNF receptor is recycled as an unbounded external receptor
	double dR_EU = dt * R_recyc_rate * pCell->custom_data[nR_IB];
	if( dR_EU > pCell->custom_data[nR_IB] )
	{ dR_EU = pCell->custom_data[nR_IB]; }

	pCell->custom_data[nR_EU] -= dR_EB; // remove newly bound receptor from R_EU 
	pCell->custom_data[nR_EB] += dR_EB; // add newly bound receptor to R_EB

	pCell->custom_data[nR_EB] -= dR_IB; // move from external bound
	pCell->custom_data[nR_IB] += dR_IB; // move to internal bound

	pCell->custom_data[nR_IB] -= dR_EU; // move from internal unbound 
	pCell->custom_data[nR_EU] += dR_EU; // move to external unbound 

	// update the TNF uptake rate 
	phenotype.secretion.uptake_rates[nTNF_external] = R_binding_rate * pCell->custom_data[nR_EU]; 

}

void update_boolean_model_input( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_TNFR" ); 
    static int nTNF_threshold = pCell->custom_data.find_variable_index( "TNFR_activation_threshold" );

    if( phenotype.death.dead == true )
	{ return; } 

    if ( pCell->custom_data[nR_EB] > pCell->custom_data[nTNF_threshold] )
	{ pCell->boolean_network.set_node_value("TNF", 1); }
	else
    { pCell->boolean_network.set_node_value("TNF", 0); }

}


void update_cell_state_model_based(Cell* pCell, Phenotype& phenotype, double dt)
{	
    static int nTNF_external = microenvironment.find_density_index( "tnf" );
    static int nTNF_export_rate = pCell->custom_data.find_variable_index( "TFN_net_production_rate" );
	static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_TNFR" ); 
    static int nTNF_threshold = pCell->custom_data.find_variable_index( "TNFR_activation_threshold" );
	static int nNFkB_activated = pCell->custom_data.find_variable_index( "NFkB_activated" );
	
	if ( pCell->boolean_network.get_node_value( "Apoptosis" ) )
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	if ( pCell->boolean_network.get_node_value( "NonACD" ) )
	{
		int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
		pCell->start_death(necrosis_model_index);
		return;
	}

	if ( pCell->boolean_network.get_node_value( "Survival" ) )
	{
		if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative )
		{ pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt); }
	}

	// If NFkB node is active produce some TNF
	double tnf_export_rate = 0;	
	if ( pCell->boolean_network.get_node_value( "NFkB" ) )
	{ 
		tnf_export_rate = pCell->custom_data[nTNF_export_rate]; 
	}
    phenotype.secretion.net_export_rates[nTNF_external] = tnf_export_rate;
}