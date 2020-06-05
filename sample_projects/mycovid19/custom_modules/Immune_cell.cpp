#include "./Immune_cell.h" 

using namespace PhysiCell; 


void Immune_cell::keep_in_bounds(double Xmin, double Xrange, double Ymin, double Yrange, double Zmin, double Zrange)
{
	static int epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	static int macrophage_type = get_cell_definition( "macrophage" ).type; 
	static int tcell_type = get_cell_definition( "CD8 Tcell" ).type; 

	if( phenotype.death.dead == false && is_out_of_domain)
	{
		is_out_of_domain = false; 
		is_active = true; 
		is_movable = true; 			
		
		std::vector<double> t_position = position; 
		t_position[0] = Xmin + Xrange * UniformRandom(); 
		t_position[1] = Ymin + Yrange * UniformRandom(); 
		
		#pragma omp critical(move_from_edge)
		{
		// std::cout << " moving cell " << this << " of type " << type_name << std::endl; 
		// std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
		
		remove_all_adhesions();
		assign_position( t_position ); 	
		update_voxel_in_container();
		
		}
	}
}    
	
	  
	  
    


