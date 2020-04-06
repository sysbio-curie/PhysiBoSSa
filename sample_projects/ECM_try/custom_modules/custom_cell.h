/* custom class for simulation that use cells that interact with extra cellular matrix

*/

#include "../core/PhysiCell_custom.h" 

#include "../BioFVM/BioFVM.h"
#include "../core/PhysiCell_phenotype.h"
#include "custom_cell_container.h"
#include "../core/PhysiCell_constants.h"
#include "../addons/PhysiBoSSa/src/boolean_network.h"
#include <math.h>
#include "../core/PhysiCell_utilities.h"

using namespace BioFVM; 
using namespace PhysiCell;

class Custom_cell : public Cell
{
private: 

Custom_cell_container * custom_container;
int freezed;

public:

    Custom_cell();
	Custom_cell_container * get_custom_container();
	void advance_bundled_phenotype_functions( double dt_ ); 
	void add_potentials(Custom_cell* other_agent);
	void custom_flag_for_removal();
	void custom_flag_for_division();
    inline bool passive() { return type == PhysiCell_constants::PASSIVE_TYPE; };
	std::vector<double> motility;
	double pintegrin;
	double pmotility;
	double padhesion;
	double nucleus_deform;
	double ecm_contact;
	double Cecm[2];
	double motility_magnitude[2];
	double Ccca_homotypic[2];
	double Ccca_heterotypic[2];
	int mmped;
	/** \brief Amount of contact with other cells */
	double cell_contact;
	/** \brief Degrade the surrounding ECM 
	*
	* @param dt time step of mechanics, to scale degradation amount
	* Currently, handle only the case of ECM as a density */
	void degrade_ecm( double dt );
	/** \brief (De)-Activate ECM degradation by the cell */
	inline void set_mmp( int activate )
	{ mmped = activate; };
    /** \brief Motility with random direction, and magnitude of motion given by customed coefficient */
	void set_3D_random_motility( double dt );
	/**
	* Motility in the polarity axis migration
	* Strength of alignement depends of the polarity parameter, as for division axis
	* Persistence defined in the polarization direction updating.
	* Polarity coefficient never reach 1 so there is some noise
	* */
	void set_3D_polarized_motility( double dt );
	void set_motility(double );

	void freezer( int frozen );
    double get_adhesion();
    bool necrotic_oxygen();
	bool has_neighbor(int);
	double adhesion(Cell* other_cell);
	double local_density(std::string field);
    /** \brief Return amount of contact with other cells */
	inline double contact_cell()
	{ return cell_contact / phenotype.geometry.radius ; };
	/** \brief Return value of adhesion strength with ECM according to integrin level */
	double integrinStrength();
	/** \brief Get the current value of heterotypic adhesion strength */
	inline double get_heterotypic_strength( double percent )
	{ return current_value( Ccca_heterotypic[0], Ccca_heterotypic[1], percent ); };
	/** \brief Get the current value of homotypic adhesion strength */
	inline double get_homotypic_strength( double percent )
	{ return current_value( Ccca_homotypic[0], Ccca_homotypic[1], percent ); };
    /** \brief Get the current value of integrin strength */
	inline double get_integrin_strength( double percent )
	{ return current_value( Cecm[0], Cecm[1], percent ); };
	/** \brief Get the current value of motility coefficient */
	inline double get_motility_amplitude( double percent )
	{ return current_value( motility_magnitude[0], motility_magnitude[1], percent ); };
	void add_cell_basement_membrane_interactions(double t, double dist);
    	/** \brief Calculate agent distance to BM if defined */
	double distance_to_membrane( double l, std::string shape);
	/** \brief Distance of agent to BA for duct geometry */
	double distance_to_membrane_duct( double l);
	/** \brief Distance of agent to BA for sphere geometry */
	double distance_to_membrane_sphere( double l);
	/** \brief Distance to membrane Sheet
	 * Basement membrane is a sheet of height 2*BM_radius 
	 * Z value is in between -BM_radius and +BM_radius
	 */
	double distance_to_membrane_sheet(double length);
	void update_cell_motion( double dt, double l, std::string shape );
	void update_velocity( double dt, double l, std::string shape );
    /** \brief Calculate repulsion and adhesion between agent and ecm at given voxel index
	 *
	 * @param index_ecm index of the ECM density in the microenv vector of densities
	 * @param index_voxel index of the current ECM voxel  */
	void add_ecm_interaction( int index_ecm, int index_voxel );
};

Custom_cell* create_custom_cell( void );  
Custom_cell* create_custom_cell( Cell_Definition& cd );  


void delete_cell( int ); 
void delete_cell( Custom_cell* ); 
void save_all_cells_to_matlab( std::string filename ); 

//function to check if a neighbor voxel contains any cell that can interact with me
bool is_neighbor_voxel(Cell* pCell, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex);