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
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
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
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	Cell_Definition* pCD = find_cell_definition( "TU" ); 
	pCD->functions.update_phenotype = TU_phenotype; 

	pCD = find_cell_definition( "Teff" ); 
	pCD->functions.update_phenotype = Teff_phenotype;


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

/*	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
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



	//place TU


	double max_TU_distance = parameters.doubles("max_TU_distance");
	Cell_Definition*pCD = find_cell_definition("TU"); 
	std::cout << "Placing cells of type " << pCD->name << "..." << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_TU") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			double r = sqrt(UniformRandom())* max_TU_distance;
			double theta = UniformRandom()*6.2831853; //*2pi*r
			position[0] = r*cos(theta); 
        	position[1] = r*sin(theta); 
        
       		 pC = create_cell( *pCD ); 
       		 pC->assign_position( position );
		}

	//place Teff
	pCD = find_cell_definition("Teff"); 
	std::cout << "Placing cells of type " << pCD->name << "..." << std::endl; 
		 for( int k=0 ; k < parameters.ints( "number_of_Teff" ); k++ )
    	{
        std::vector<double> position = UniformOnUnitCircle(); 
        position *= parameters.doubles("max_Teff_distance"); 
        
        pC = create_cell( *pCD ); 
        pC->assign_position( position );
   		 }


	//place Texh
	pCD = find_cell_definition("Texh"); 
	std::cout << "Placing cells of type " << pCD->name << "..." << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_Texh") ; n++ )
		{
        std::vector<double> position = UniformOnUnitCircle(); 
        position *= parameters.doubles("max_Teff_distance"); 
        
        pC = create_cell( *pCD ); 
        pC->assign_position( position );
   		 }

	//place M1 Mph
	double max_Mph_distance = parameters.doubles("max_Mph_distance");
	pCD = find_cell_definition("M1 Mph"); 
	std::cout << "Placing cells of type " << pCD->name << "..." << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_M1") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			double r = sqrt(UniformRandom())* max_Mph_distance;
			double theta = UniformRandom()*6.2831853; //*2pi*r
			position[0] = r*cos(theta); 
        	position[1] = r*sin(theta); 
        
       		 pC = create_cell( *pCD ); 
       		 pC->assign_position( position );
		}

	//place M2 Mph
	pCD = find_cell_definition("M2 Mph"); 
	std::cout << "Placing cells of type " << pCD->name << "..." << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_M2") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			double r = sqrt(UniformRandom())* max_Mph_distance;
			double theta = UniformRandom()*6.2831853; //*2pi*r*
			position[0] = r*cos(theta); 
        	position[1] = r*sin(theta); 
        
       		 pC = create_cell( *pCD ); 
       		 pC->assign_position( position );
		}
*/

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

//phenotype of Tumor cells

void TU_phenotype(Cell*pCell , Phenotype& phenotype, double dt)
{
	static Cell_Definition* pCD = find_cell_definition (pCell -> type_name); 


	static int nApoptosis = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model); 

	double damage_signal = pCell -> state.damage; 
	double base_death_rate = pCD->phenotype.death.rates [nApoptosis]; 
	double max_death_rate = 0.05; 

	static double damage_halfmax = 36.0; 
	double hill_power=2.0; 

	double hill = Hill_response_function(damage_signal, damage_halfmax, hill_power); 
	phenotype.death.rates[nApoptosis] = base_death_rate + (max_death_rate)*hill;

	return;
}

//phenotype of Teff

void Teff_phenotype(Cell*pCell , Phenotype& phenotype, double dt)
{
//AI induces transformation to Texh 
	double max_transform_rate = 0.0083; 

	double AI = get_single_signal(pCell, "AI");
	double AI_halfmax = 2;
	double AI_hill_power=5;

	double hill_transform_rate = Hill_response_function(AI, AI_halfmax, AI_hill_power); 
	double transform_rate = get_single_behavior (pCell, "transform to Texh") + max_transform_rate*hill_transform_rate;
	set_single_behavior(pCell, "transform to Texh", transform_rate);

//if attached to TU and TU is dead, then detach and transform to Texh

	if (pCell -> state.number_of_attached_cells () > 0.5)
	{	
		
		for( int k=0; k < pCell->state.attached_cells.size() ; k++ )
		{
			Cell* pTemp = pCell -> state.attached_cells[k]; 

			//if attached to TU (defined with custom receptor = 1), then stop migration and proliferation 
		
			double receptor_attached = get_single_signal(pTemp, "custom:receptor"); 
			if (receptor_attached > 0.5)
			{
				set_single_behavior(pCell, "migration speed", 0.0);
				set_single_behavior(pCell, "cycle entry", 0.0);

				//if TU is killed, then transform to Texh and detach
				if (pTemp->phenotype.death.dead == true)
				{
				set_single_behavior(pCell, "transform to Texh", 1.0);
				pCell -> remove_all_attached_cells();
				}
			}
		}

	}	

//if unattached and close to TU, attach to TU and attack dependent on PI

	static double elastic_coefficient = parameters.doubles ("elastic_coefficient"); 

	set_single_behavior (pCell, "cell-cell adhesion elastic constant", elastic_coefficient); 

	double attack_rate = get_single_behavior(pCell, "attack TU"); 

	if (pCell-> state.number_of_attached_cells() == 0)
	{
		std::vector<Cell*> nearby = pCell -> cells_in_my_container();
		for (int i=0; i < nearby.size(); i++)
		{
			double receptor = get_single_signal (nearby[i], "custom:receptor");

			if (receptor > 0.5 && nearby[i]->phenotype.death.dead == false )
			{
				attach_cells (pCell, nearby[i]); 
				double max_attack_rate = 100.0; 
				double PI = get_single_signal(pCell, "PI"); 
				double PI_halfmax = 0.5; 
				double PI_hill_power = 5; 
				double hill_attack_rate = Hill_response_function (PI, PI_halfmax, PI_hill_power); 
				double attack_rate = get_single_behavior (pCell, "attack TU") + max_attack_rate*hill_attack_rate; 
				set_single_behavior(pCell, "attack TU", attack_rate); 
				
			}
		}
	} 

	return;

}
//coloring function

std::vector<std::string> custom_coloring_function( Cell* pCell )
{
	// start with color-by-number (as above)
	std::vector<std::string> output = paint_by_number_cell_coloring(pCell);
	// dead cancer cells: brown 
	
	if( pCell->type_name == "TU") 
	{ bool dead = (bool) get_single_signal( pCell, "dead"); 
	if( dead )
	{
	output[0] = "rgb(111,78,55)";
	}
	return output;
	}
}

