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

#include "./custom.h"

// declare cell definitions here 

Cell_Definition fibroblast; 
Cell_Definition organoid_cell;


void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "default cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// only needed for a 2-D simulation: 
	
	/*
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	*/
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int Start_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model);
	int End_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model );
	

	// cell_defaults no necrosis and apoptosis
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.cycle.data.transition_rate(Start_index,End_index) = 0.0;


	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0; 
	
	// add custom data here, if any 
	
	
	
	// --------- Defining Organoid Cells -------- //
	organoid_cell = cell_defaults; 
	organoid_cell.type = 1; 
	organoid_cell.name = "organoid cell"; 
	
	// make sure the new cell type has its own reference phenotyhpe
	
	organoid_cell.parameters.pReference_live_phenotype = &( organoid_cell.phenotype ); 
	
	organoid_cell.phenotype.motility.is_motile = false; 
	// organoid_cell.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; // 15 minutes
	// organoid_cell.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25; // 0.25 micron/minute 
	// organoid_cell.phenotype.motility.migration_bias = 0.0;// completely random 
	
    organoid_cell.phenotype.molecular.sync_to_microenvironment( &microenvironment );
    organoid_cell.functions.update_phenotype = tumor_energy_update_function;
    
    
	// Set cell-cell adhesion to 5% of other cells 
	organoid_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles( "organoid_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	organoid_cell.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "organoid_cell_apoptosis_rate" ); // 0.0; 
	
	// Setting Proliferation
	// 
	organoid_cell.phenotype.cycle.data.transition_rate(Start_index,End_index) = parameters.doubles( "organoid_cell_relative_cycle_entry_rate" ); // 0.1; 
		
		
		
	//------------ Defining Fibroblast ------------//
	
	fibroblast = cell_defaults; 
	fibroblast.type = 2; 
	fibroblast.name = "fibroblast"; 
	
	// make sure the new cell type has its own reference phenotype
	
	fibroblast.parameters.pReference_live_phenotype = &( fibroblast.phenotype ); 
	// Set cell-cell adhesion to 5% of other cells 
	//fibro_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "fibroblast_relative_adhesion" );
	
	// Set apoptosis to zero 
	fibroblast.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "fibroblast_apoptosis_rate" ); // 0.0; 
	fibroblast.phenotype.cycle.data.transition_rate(Start_index,End_index) = 0.0;
	fibroblast.phenotype.cycle.data.transition_rate(End_index,Start_index) = 0.0;
	
	fibroblast.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0;
	fibroblast.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0; 
	fibroblast.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0; 
	
	//fibroblast.phenotype.secretion.uptake_rates[glucose_substrate_index] = 0.0; 
	//fibroblast.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0.0; 
	//fibroblast.phenotype.secretion.saturation_densities[glucose_substrate_index] = 0.0; 
	
	//fibroblast.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; 
	//fibroblast.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.0; 
	//fibroblast.phenotype.secretion.saturation_densities[lactate_substrate_index] = 0.0; 
	
	// ---- END -- Fibroblast Cell Definitions -- END ---- //	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/*	
	default_microenvironment_options.X_range = {-500, 500}; 
	default_microenvironment_options.Y_range = {-500, 500}; 
	default_microenvironment_options.Z_range = {-500, 500}; 
*/	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	
	
	
/*
	// all this is in XML as of August 2019 (1.6.0)
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pCell;

    if (parameters.bools("fibroblast_seeding"))
        {
            std::cout << "creating fibroblasts" << std::endl;
		
            for (int i= -2666; i<2666; i+=16.82)
            {
                for (int j= -2666; j<2666; j+=16.82)
                {			
                    pCell = create_cell(fibroblast);
                    pCell->assign_position(i,-500,j);	
                }
            } 	
        }
	
	if (parameters.bools("organoid_cell_seeding"))
		{
			std::cout << "creating organoid cells" << std::endl;
			SeedRandom( parameters.ints("random_seed") ); // or specify a seed here
			for (int i = 0; i < parameters.ints("organoid_cell_number"); i++) // seeding number of organoid cells specified in PhysiCell_settings.xml
			{
				int x_pos, y_pos, z_pos; // correspond to positions on x, y, and z axes
				x_pos = (rand() % 5333) - 2666; 
				y_pos = (rand() % 961) - 480;
				z_pos = (rand() % 5333) - 2666;
				pCell = create_cell(organoid_cell);
				pCell->assign_position(x_pos, y_pos, z_pos);
			}
		}	

	return; 
}

void tumor_energy_update_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	static int Start_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model);
	static int End_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model );

    double tr = phenotype.cycle.data.transition_rate( Start_index,End_index );
    //std::cout << tr << std::endl;
    
	return;
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// if the cell is cell and not dead, paint it black 
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 1 )
	{
		 output[0] = "red"; 
		 output[2] = "red"; 	
	}
	
	return output; 
}
