#ifndef __Custom_cell_h__
#define __Custom_cell_h__

/* custom class for simulation that use cells that interact with extra cellular matrix
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include "../core/PhysiCell_cell.h" 
#include "../core/PhysiCell_constants.h"
#include "custom_main.h"
using namespace PhysiCell;

class Custom_cell : public Cell {

private:
	int freezed;	
	double ecmrad;
	
protected:
	
public:
	inline bool passive() { 
		return (
			custom_data.find_variable_index("passive_cell") >= 0 
			&& ((int) custom_data["passive_cell"]) == 1
		); 
	};
	
    std::vector<double> motility;
	double pintegrin;
	double pmotility;
	double padhesion;
	double nucleus_deform;
	double ecm_contact;
	int mmped;
	Custom_cell();


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

    /** \brief Return amount of contact with other cells */
	inline double contact_cell()
	{ return cell_contact / phenotype.geometry.radius ; };
	
	/** \brief Return value of adhesion strength with ECM according to integrin level */
	double integrinStrength();
	
	/** \brief Get the current value of heterotypic adhesion strength */
	inline double get_heterotypic_strength( double percent )
	{ return current_value( PhysiCell::parameters.doubles("heterotypic_adhesion_min"), PhysiCell::parameters.doubles("heterotypic_adhesion_max"), percent ); };
	
	/** \brief Get the current value of homotypic adhesion strength */
	inline double get_homotypic_strength( double percent )
	{ return current_value( PhysiCell::parameters.doubles("homotypic_adhesion_min"), PhysiCell::parameters.doubles("homotypic_adhesion_max"), percent ); };

    /** \brief Get the current value of integrin strength */
	inline double get_integrin_strength( double percent )
	{ return current_value( PhysiCell::parameters.doubles("ecm_adhesion_min"), PhysiCell::parameters.doubles("ecm_adhesion_max"), percent ); };
	
	/** \brief Get the current value of motility coefficient */
	inline double get_motility_amplitude( double percent )
	{ return current_value(PhysiCell::parameters.doubles("motility_amplitude_min"), PhysiCell::parameters.doubles("motility_amplitude_max"), percent ); };
    
    bool has_neighbor(int);
	double adhesion(Cell* other_cell);
	double get_adhesion();

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

    /** \brief Calculate repulsion and adhesion between agent and ecm at given voxel index
	 *
	 * @param index_ecm index of the ECM density in the microenv vector of densities
	 * @param index_voxel index of the current ECM voxel  */
	void add_ecm_interaction( int index_ecm, int index_voxel );

	static Cell* create_custom_cell();

	static void check_passive(Cell* cell, Phenotype& phenotype, double dt);

	static void custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt);
	static double custom_repulsion_function(Cell* pCell, Cell* otherCell);
	static double custom_adhesion_function(Cell* pCell, Cell* otherCell, double distance);
	static bool wait_for_nucleus_growth(Cell* cell, Phenotype& phenotype, double dt);
	static bool waiting_to_remove(Cell* cell, Phenotype& phenotype, double dt);
	/** \brief Change the current value of motility percent coeff, increase or decrease according to up value */
	inline void evolve_motility_coef( int up, double dt )
	{ evolve_coef( up, &pmotility, dt); };
	/** \brief Set the value of freezed */
	inline void freezing( int frozen )
	{ freezed = frozen; };
	/** \brief Return amount of contact with ECM */
	inline double contact_ecm()
	{ return ecm_contact / phenotype.geometry.radius ; };
	void set_initial_volume(Cell_Definition cell_def, float radius);
};

#endif