

#include "./tnf_response_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information tnf_dynamics_info;

void tnf_dynamics_model_setup()
{
    tnf_dynamics_info.name = "Tumor Necrotic Factor model dynamics"; 
	tnf_dynamics_info.version = "0.1.0";
	
    tnf_dynamics_info.main_function= tnf_dynamics_model; 

	// what custom data do I need? 
    tnf_dynamics_info.cell_variables.push_back( "unbound external TNFR" );
	tnf_dynamics_info.cell_variables.push_back( "bound external TNFR" );
	tnf_dynamics_info.cell_variables.push_back( "bound internal TNFR" );
    
    tnf_dynamics_info.cell_variables.push_back( "TNFR binding rate" ); 
	tnf_dynamics_info.cell_variables.push_back( "TNFR endocytosis rate" );
    tnf_dynamics_info.cell_variables.push_back( "TNFR recycling rate" );
    tnf_dynamics_info.cell_variables.push_back( "TNFR activation threshold" );
    tnf_dynamics_info.cell_variables.push_back( "NFkB activated" );
	tnf_dynamics_info.cell_variables.push_back( "TFN net production rate" );

	tnf_dynamics_info.register_model();
}


void tnf_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nTNF_external = microenvironment.find_density_index( "tnf" );
    static int nR_EU = pCell->custom_data.find_variable_index( "unbound external TNFR" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound external TNFR" );
	static int nR_IB = pCell->custom_data.find_variable_index( "bound internal TNFR" );

    static int nR_bind = pCell->custom_data.find_variable_index( "TNFR binding rate" ); 

    static int nR_endo = pCell->custom_data.find_variable_index( "TNFR endocytosis rate" ); 
	static int nR_recycle = pCell->custom_data.find_variable_index( "TNFR recycling rate" ); 

	if( phenotype.death.dead == true )
	{ return; } 
		
    // internalized TNF tells us how many have recently bound to receptors
	// TNF is internalized at:
	// phenotype.secretion.uptake_rates[nTNF_external] = 
	// 					pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU];
	
	// The internalization is only used to track the TNF
	// The following part of the code takes care of correcly managed
	double newly_bound = phenotype.molecular.internalized_total_substrates[nTNF_external]; 
	
	// if it tried to bind to more virus than there are receptors, compensate 
	double excess_binding = newly_bound - pCell->custom_data[nR_EU]; 
	if( excess_binding > 0.0 )
	{
		// don't bring in more virus than there are receptors 
		newly_bound = pCell->custom_data[nR_EU]; 
		// dump any excess back into the microenvironment
		static double to_density = 1.0 / microenvironment.mesh.dV; 
		// this needs omp critical because 2 cells writing to 1 voxel is not thread safe 
		#pragma omp critical 
		{
			pCell->nearest_density_vector()[nTNF_external] += excess_binding * to_density; 
		}
	}
	phenotype.molecular.internalized_total_substrates[nTNF_external] = 0.0; 
	
	// add newly bound receptor to R_EB
	pCell->custom_data[nR_EB] += newly_bound; 
	
	// remove newly bound receptor from R_EU 
	pCell->custom_data[nR_EU] -= newly_bound; 
	
	// endocytosis 
	double dR_IB = dt * pCell->custom_data[nR_endo] * pCell->custom_data[nR_EB];
	
    if( dR_IB > pCell->custom_data[nR_EB] )
	{ dR_IB = pCell->custom_data[nR_EB]; }
	
    pCell->custom_data[nR_EB] -= dR_IB; // move from external bound
	pCell->custom_data[nR_IB] += dR_IB; // move to internal bound
	
	// TNF release from receptor 
	
	double dR_IU = dt * pCell->custom_data[nR_recycle] * pCell->custom_data[nR_IB];
	
    if( dR_IU > pCell->custom_data[nR_IB] )
	{ dR_IU = pCell->custom_data[nR_IB]; }
	
    // The internalized bounded TNFR release the TNF
    // The TNF is instantaneously degraded by the cell
    // The TNF receptor is recycled as an 
    pCell->custom_data[nR_IB] -= dR_IU; // move from internal bound 
    pCell->custom_data[nR_EU] += dR_IU; // move to internal unbound 
	
	// receptor recycling 
	double dR_EU = dt*pCell->custom_data[nR_recycle]*pCell->custom_data[nR_IB];
	if( dR_EU > pCell->custom_data[nR_IB] )
	{ dR_EU = pCell->custom_data[nR_IB]; }
	pCell->custom_data[nR_IB] -= dR_EU; // move from internal unbound 
	pCell->custom_data[nR_EU] += dR_EU; // move to external unbound 
	
	// update the virion uptake rate 
	phenotype.secretion.uptake_rates[nTNF_external] = pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU]; 

}

void update_boolean_model_input( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int nR_EB = pCell->custom_data.find_variable_index( "bound external TNFR" ); 
    static int nTNF_threshold = pCell->custom_data.find_variable_index( "TNFR activation threshold" );

    if( phenotype.death.dead == true )
	{ return; } 

    if ( pCell->custom_data[nR_EB] > pCell->custom_data[nTNF_threshold] )
    { 
		std::cout << "TNF Activated! "  << std::endl;
		pCell->boolean_network.set_node_value("TNF", 1); }
	else
    { pCell->boolean_network.set_node_value("TNF", 0); }

}

void update_custom_variables( Cell* pCell )
{
    static int tnf_external = microenvironment.find_density_index( "tnf" ); 
	static int tnf_internal = pCell->custom_data.find_variable_index( "tnf" );
    static int index_tnf_node = pCell->custom_data.find_variable_index("tnf_node");
	static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
		
	pCell->custom_data.variables.at(index_tnf_node).value = pCell->boolean_network.get_node_value("TNF");
	pCell->custom_data.variables.at(index_fadd_node).value = pCell->boolean_network.get_node_value("FADD");
}

void update_cell_state_model_based(Cell* pCell, Phenotype& phenotype, double dt)
{


    static int nTNF_external = microenvironment.find_density_index( "tnf" );
    static int nTNF_export_rate = pCell->custom_data.find_variable_index( "TFN net production rate" );
    double tnf_export_rate = 0;

	static int nR_EB = pCell->custom_data.find_variable_index( "bound external TNFR" ); 
    static int nTNF_threshold = pCell->custom_data.find_variable_index( "TNFR activation threshold" );
	bool tnf_active = pCell->custom_data[nR_EB] > pCell->custom_data[nTNF_threshold];

	if ( pCell->boolean_network.get_node_value( "Apoptosis" ) )
	{
		
		std::cout << "Apoptosis! " << tnf_active << std::endl;
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		return;
	}

	if ( pCell->boolean_network.get_node_value( "NonACD" ) )
	{
		std::cout << "NonACD! " << tnf_active << std::endl;
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
	if ( pCell->boolean_network.get_node_value( "NFkB" ) )
	{ 
		std::cout << "NFkB! " << tnf_active << std::endl;
		tnf_export_rate = pCell->custom_data[nTNF_export_rate]; 
	}
    
    phenotype.secretion.net_export_rates[nTNF_external] = tnf_export_rate;
}