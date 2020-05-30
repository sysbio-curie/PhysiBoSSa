#include "maboss_intracellular.h"

MaBoSSIntracellular::MaBoSSIntracellular() : Intracellular()
{
	type = "maboss";
	initial_values.clear();
	mutations.clear();
	parameters.clear();
}

MaBoSSIntracellular::MaBoSSIntracellular(pugi::xml_node& node)
{
	// MaBoSSIntracellular* pIntra = new MaBoSSIntracellular();
	type = "maboss";
	initialize_intracellular_from_pugixml(node);

}

MaBoSSIntracellular::MaBoSSIntracellular(MaBoSSIntracellular* copy) 
{
		type = copy->type;
		bnd_filename = copy->bnd_filename;
		cfg_filename = copy->cfg_filename;
		time_step = copy->time_step;
		std::map<std::string, double> initial_values(copy->initial_values);
		std::map<std::string, double> mutations(copy->mutations);
		std::map<std::string, double> parameters(copy->parameters);
		
		if (copy->network.has_init()) {
			network.initialize_boolean_network(
				copy->bnd_filename, copy->cfg_filename,
				copy->initial_values, copy->mutations, copy->parameters
			);
			network.set_time_step(copy->time_step);
			network.set_discrete_time(copy->discrete_time, copy->time_tick);
		}
		
}

void MaBoSSIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{
	pugi::xml_node node_bnd = node.child( "bnd_filename" );
	if ( node_bnd )
	{ bnd_filename = PhysiCell::xml_get_my_string_value (node_bnd); }
	
	pugi::xml_node node_cfg = node.child( "cfg_filename" );
	if ( node_cfg )
	{ cfg_filename = PhysiCell::xml_get_my_string_value (node_cfg); }
	
	pugi::xml_node node_init_values = node.child( "initial_values" );
	if( node_init_values )
	{
		pugi::xml_node node_init_value = node_init_values.child( "initial_value" );
		while( node_init_value )
		{
			std::string node_name = node_init_value.attribute( "node" ).value(); 
			double node_value = PhysiCell::xml_get_my_double_value( node_init_value );
			
			initial_values[node_name] = node_value;
			
			node_init_value = node_init_value.next_sibling( "initial_value" ); 
		}
	}
	
	pugi::xml_node node_mutations = node.child( "mutations" );
	if( node_mutations )
	{
		pugi::xml_node node_mutation = node_mutations.child( "mutation" );
		while( node_mutation )
		{
			std::string node_name = node_mutation.attribute( "node" ).value(); 
			double node_value = PhysiCell::xml_get_my_double_value( node_mutation );
			
			mutations[node_name] = node_value;
			
			node_mutation = node_mutation.next_sibling( "mutation" ); 
		}
	}
	
	pugi::xml_node node_parameters = node.child( "parameters" );
	if( node_parameters )
	{
		pugi::xml_node node_parameter = node_parameters.child( "parameter" );
		while( node_parameter )
		{
			std::string param_name = node_parameter.attribute( "name" ).value(); 
			double param_value = PhysiCell::xml_get_my_double_value( node_parameter );
			
			parameters[param_name] = param_value;
			
			node_parameter = node_parameter.next_sibling( "parameter" ); 
		}
	}
	
	network.initialize_boolean_network(
		bnd_filename, 
		cfg_filename,
		initial_values,
		mutations,
		parameters
	);
	
	pugi::xml_node node_timestep = node.child( "time_step" ); 
	if( node_timestep )
	{ 
		time_step = PhysiCell::xml_get_my_double_value( node_timestep );
		network.set_time_step(time_step);
	}
	
	pugi::xml_node node_discretetime = node.child( "discrete_time" ); 
	pugi::xml_node node_timetick = node.child( "time_tick" ); 

	if( node_discretetime && node_timetick )
	{ 
		discrete_time = PhysiCell::xml_get_my_bool_value( node_discretetime );		
		time_tick = PhysiCell::xml_get_my_double_value( node_timetick );
		network.set_discrete_time(discrete_time, time_tick);
	}

	
	
}

MaBoSSIntracellular* getMaBoSSModel(PhysiCell::Phenotype& phenotype) {
	return static_cast<MaBoSSIntracellular*>(phenotype.intracellular);
}

