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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

void create_cell_types( void )
{
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	Cell_Definition* pCD = find_cell_definition( "Tcell" ); 
	pCD->functions.update_phenotype = Teff_phenotype; 

	pCD = find_cell_definition( "Texh" );
	pCD->functions.update_phenotype = Texh_phenotype;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 

	// initialize BioFVM 
	
	initialize_microenvironment(); 	

	double cyt_a_conc = parameters.doubles("p_cyt_a"); //get concentration of cyt_a from user_parameters

	std::vector<double> int_cond = {0, 0, 0, 0, cyt_a_conc}; //define values for initial conditions {T_chem, Mph_chem, T_cyt, debris, cyt_a}
	for (int i = 0; i < microenvironment.number_of_voxels(); i++)
	{
		microenvironment.add_dirichlet_node(i, int_cond);
	}
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

//Phenotype of Tcell (effector Tcells)
void Teff_phenotype(Cell*pCell , Phenotype& phenotype, double dt)
{
//Tcell exhaustion (simulated by transformation of Tcell to Texh)
//contact to Mph_TU or Mph_IM induces transformation to Texh, if Mph_TU_AI = 1 or Mph_IM_AI=1 (#1)
//attack time with TU induces transformation to Texh, if CP = 1 on Tcells (#2)
//intracellular amount of of cytokine (IL6) induces transformation to Texh (#3)

	double Mph_TU_AI = parameters.doubles("Mph_TU_AI"); 
	double Mph_IM_AI = parameters.doubles("Mph_IM_AI"); 
	double CP = parameters.doubles("CP"); 

	double contact_Mph = get_single_signal(pCell, "contact with Mph_TU")*Mph_TU_AI + get_single_signal(pCell, "contact with Mph_IM")*Mph_IM_AI; //#1 in vector, contact with Mph_TU or Mph_IM affects Tcell exhaustion, depndent on Mph_TU_AI or Mph_IM_AI expression (simulates secretion of anti-inflammatory factors by Mph)
	double attack_TU = (get_single_signal(pCell, "damage delivered")/60)*CP; //#2 in vector, attacking tumor cells affects Tcell exhaustion, dependent on CP expression (simulates checkpoint expression on Tcells), one attack duration = 60min, attack_damage_rate = 1/min, so one attack delivers 60 damage, /60 for one attack
	double cyt_a = get_single_signal(pCell, "intracellular cyt_a"); //#3 in vetor

	std::vector<double> exhaustion_parameter = {contact_Mph, attack_TU, cyt_a}; //exhaustion parameter combined value of contact with Mph/EC and attack with tumor cells

	double max_transform_rate = parameters.doubles("Tcell_transform_max"); 
	double cyt_a_halfmax = parameters.doubles("cyt_a_halfmax"); //halfmax for cyt_a

	std::vector<double> transform_halfmax = {5,5, cyt_a_halfmax}; //halfmaxes for contact with Mph, attack with TU and IL6
	std::vector<double> transform_hillpower = {20,20,4}; //hill powers for contact with Mph, attack with TU and IL6

	double hill_transform_rate = multivariate_Hill_response_function(exhaustion_parameter, transform_halfmax, transform_hillpower); //multivariate Hill function for exhaustion parameter
	double transform_rate = get_single_behavior (pCell, "transform to Texh") + max_transform_rate*hill_transform_rate;

	set_single_behavior(pCell, "transform to Texh", transform_rate); 

//Tcell arrest
//if unattached and close to Mph, attach to Mph and stop migration and proliferation, if chem_r on Mph is expressed (chem_r=1)

	if (pCell-> state.number_of_attached_cells() == 0)
	{
		std::vector<Cell*> nearby = pCell -> cells_in_my_container();
		for (int i=0; i < nearby.size(); i++)
		{
			double receptor = get_single_signal (nearby[i], "custom:chem_r");

			if (receptor == 1)
			{
				attach_cells (pCell, nearby[i]); 
				set_single_behavior(pCell, "migration speed", 0.0);
				set_single_behavior(pCell, "cycle entry", 0.0);
				/*if (exh == 1)
				{
				pCell -> remove_all_attached_cells();
				}*/
			}
		}
	}
	return;	
}


//Phenotype of Texh
void Texh_phenotype(Cell*pCell , Phenotype& phenotype, double dt)
{
//Tcell arrest
//if unattached and close to Mph, attach to Mph and stop migration and proliferation, if chem_r on Mph is expressed (chem_r=1)

	if (pCell-> state.number_of_attached_cells() == 0)
	{
		std::vector<Cell*> nearby = pCell -> cells_in_my_container();
		for (int i=0; i < nearby.size(); i++)
		{
			double receptor = get_single_signal (nearby[i], "custom:chem_r");

			if (receptor == 1)
			{
				attach_cells (pCell, nearby[i]); 
				set_single_behavior(pCell, "migration speed", 0.2);
				set_single_behavior(pCell, "cycle entry", 0.0);
				/*if (exh == 1)
				{
				pCell -> remove_all_attached_cells();
				}*/
			}

		}
	} 
	return;
}

/*
//coloring function

std::vector<std::string> custom_coloring_function(Cell* pCell )
{
	// start with color-by-number (as above)
	std::vector<std::string> output = paint_by_number_cell_coloring(pCell);
	// dead cancer cells: brown 
	bool apoptotic = (bool) get_single_signal(pCell, "apoptotic");

	if( pCell->type_name == "TU" && apoptotic > 0.5) 
	{ 
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}
	return output;
}
*/	


