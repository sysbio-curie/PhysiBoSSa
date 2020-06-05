#include "./Attachable_cell.h" 

using namespace PhysiCell; 
    
void Attachable_cell::attach( Cell* pCell )
{
	#pragma omp critical(attach)
	{
		bool already_attached = false; 
		for( int i=0 ; i < state.neighbors.size() ; i++ )
		{
			if( state.neighbors[i] == pCell )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ state.neighbors.push_back( pCell ); }
		
		already_attached = false; 
		for( int i=0 ; i < pCell->state.neighbors.size() ; i++ )
		{
			if( pCell->state.neighbors[i] == this )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ pCell->state.neighbors.push_back( this ); }
	}

	return; 
}

void Attachable_cell::detach( Cell* pCell )
{
	#pragma omp critical(detach)
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( state.neighbors[i] == pCell )
			{
				int n = state.neighbors.size(); 
				// copy last entry to current position 
				state.neighbors[i] = state.neighbors[n-1]; 
				// shrink by one 
				state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell->state.neighbors[i] == this )
			{
				int n = pCell->state.neighbors.size(); 
				// copy last entry to current position 
				pCell->state.neighbors[i] = pCell->state.neighbors[n-1]; 
				// shrink by one 
				pCell->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	
	return; 
}

	  
void Attachable_cell::remove_all_adhesions()
{
	// detach all attached cells 
	for( int n = 0; n < state.neighbors.size() ; n++ )
	{ detach( state.neighbors[n] ); }		
	
	return; 
}



