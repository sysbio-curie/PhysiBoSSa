/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#ifndef __PhysiCell_cell_h__
#define __PhysiCell_cell_h__

#include "./PhysiCell_custom.h" 

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_phenotype.h"
#include "./PhysiCell_cell_container.h"
#include "./PhysiCell_constants.h"
#include "../addons/PhysiBoSSa/src/boolean_network.h"
#include <math.h>
#include "../core/PhysiCell_utilities.h"
using namespace BioFVM; 

namespace PhysiCell{
class Cell_Container;

class Cell_Parameters
{
 private:
 public:
	// oxygen values (in mmHg) for critical phenotype changes
	double o2_hypoxic_threshold; // value at which hypoxic signaling starts
	double o2_hypoxic_response; // value at which omics changes are observed 
	double o2_hypoxic_saturation; // value at which hypoxic signalign saturates 
	// o2_hypoxic_saturation < o2_hypoxic_threshold
	
	double o2_proliferation_saturation; // value at which extra o2 does not increase proliferation
	double o2_proliferation_threshold; // value at which o2 is sufficient for proliferation

	double o2_reference; // physioxic reference value, in the linked reference Phenotype
	// o2_proliferation_threshold < o2_reference < o2_proliferation_saturation; 
	
	double o2_necrosis_threshold; // value at which cells start experiencing necrotic death 
	double o2_necrosis_max; // value at which necrosis reaches its maximum rate 
	// o2_necrosis_max < o2_necrosis_threshold
	
	Phenotype* pReference_live_phenotype; // reference live phenotype (typically physioxic) 
	// Phenotype* pReference_necrotic_phenotype; // reference live phenotype (typically physioxic) 

	// necrosis parameters (may evenually be moved into a reference necrotic phenotype 
	double max_necrosis_rate; // deprecate
	int necrosis_type; // deprecate 
	
	
	Cell_Parameters(); 
}; 

class Cell_Definition
{
 private:
 public: 
	int type; 
	std::string name; 
 
	Microenvironment* pMicroenvironment; 
	
	Cell_Parameters parameters; 
	Custom_Cell_Data custom_data; 
	Cell_Functions functions; 
	Phenotype phenotype; 

	Cell_Definition();  // done 
	Cell_Definition( Cell_Definition& cd ); // copy constructor 
	Cell_Definition& operator=( const Cell_Definition& cd ); // copy assignment 
};

extern Cell_Definition cell_defaults; 

class Cell_State
{
 public:
	std::vector<Cell*> neighbors; // not currently tracked! 
	std::vector<double> orientation;
	
	double simple_pressure; 
	
	Cell_State(); 
};

class Cell : public Basic_Agent 
{
 private: 
	Cell_Container * container;
	int current_mechanics_voxel_index;
	int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index
	int freezed;	
 public:
	std::string type_name; 
 
	Custom_Cell_Data custom_data;
	Cell_Parameters parameters;
	Cell_Functions functions; 
	inline bool passive() { return type == PhysiCell_constants::PASSIVE_TYPE; };
	Cell_State state; 
	Phenotype phenotype; 
	
	BooleanNetwork boolean_network;
	
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
	void update_motility_vector( double dt_ );
	void advance_bundled_phenotype_functions( double dt_ ); 
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
	void add_potentials(Cell*);       // Add repulsive and adhesive forces.
	void set_previous_velocity(double xV, double yV, double zV);
	int get_current_mechanics_voxel_index();
	void turn_off_reactions(double); 		  // Turn off all the reactions of the cell
	void freezer( int frozen );
	bool is_out_of_domain;
	bool is_movable;
	double get_adhesion();
	void flag_for_division( void ); // done 
	void flag_for_removal( void ); // done 
	
	void start_death( int death_model_index ); 
	void lyse_cell( void ); 

	Cell* divide( void );
	void die( void );
	void step(double dt);
	Cell();
	
	bool assign_position(std::vector<double> new_position);
	bool assign_position(double, double, double);
	void set_total_volume(double);
	bool necrotic_oxygen();
	bool has_neighbor(int);
	double adhesion(Cell* other_cell);
	double local_density(std::string field);
	double& get_total_volume(void); // NEW
	// mechanics 
	void update_position( double dt ); //
	std::vector<double> displacement; // this should be moved to state, or made private  

	
	void assign_orientation();  // if set_orientaion is defined, uses it to assign the orientation
								// otherwise, it assigns a random orientation to the cell.
	
	void copy_function_pointers(Cell*);
	
	void update_voxel_in_container(void);
	void copy_data(Cell *);
	
	void ingest_cell( Cell* pCell_to_eat ); // for use in predation, e.g., immune cells 
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
	// I want to eventually deprecate this, by ensuring that 
	// critical BioFVM and PhysiCell data elements are synced when they are needed 
	/** \brief Get the current value of integrin strength */
	inline double get_integrin_strength( double percent )
	{ return current_value( Cecm[0], Cecm[1], percent ); };
	/** \brief Get the current value of motility coefficient */
	inline double get_motility_amplitude( double percent )
	{ return current_value( motility_magnitude[0], motility_magnitude[1], percent ); };
	void add_cell_basement_membrane_interactions(double t, double dist);
	void set_phenotype( Phenotype& phenotype ); // no longer needed?
	void update_radius();
	Cell_Container * get_container();
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
	std::vector<Cell*>& cells_in_my_container( void ); 
			/** \brief Calculate repulsion and adhesion between agent and ecm at given voxel index
	 *
	 * @param index_ecm index of the ECM density in the microenv vector of densities
	 * @param index_voxel index of the current ECM voxel  */
	void add_ecm_interaction( int index_ecm, int index_voxel );

	void convert_to_cell_definition( Cell_Definition& cd ); 
};

Cell* create_cell( void );  
Cell* create_cell( Cell_Definition& cd );  


void delete_cell( int ); 
void delete_cell( Cell* ); 
void save_all_cells_to_matlab( std::string filename ); 

//function to check if a neighbor voxel contains any cell that can interact with me
bool is_neighbor_voxel(Cell* pCell, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex);  

};

#endif