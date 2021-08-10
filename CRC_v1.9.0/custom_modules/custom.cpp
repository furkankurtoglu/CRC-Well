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

int oxygen_i, glucose_i; 
int energy_vi; 

// These are for C
//#define STATIC_RRC
// #include "rrc_api.h"
// #include "rrc_types.h"
// #include "rrc_utilities.h"
// extern "C" rrc::RRHandle createRRInstance();

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// #include <vector>
#include <string>


void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();

  	initialize_cell_definitions_from_pugixml();   //rwh

    //rwh
    cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	// cell_defaults.functions.update_migration_bias = weighted_motility_function; 
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

    build_cell_definitions_maps();
    display_cell_definitions(std::cout);
    /*fibroblast.functions.update_phenotype = tumor_energy_update_function_fibroblast; */
    /*KRAS_positive.functions.update_phenotype = tumor_energy_update_function_KRAS_positive;*/
    //KRAS_negative.functions.update_phenotype = tumor_energy_update_function_KRAS_negative;
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
	
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}


void setup_1D_microenvironment( void )
{
    Microenvironment coarse_well;
    std::cout << "SPOT 1" <<std::endl;
    
    coarse_well.name = "coarse well";
    coarse_well.spatial_units = "micron";
    coarse_well.mesh.units = "micron";
    coarse_well.time_units = "min";

    coarse_well.set_density( 0 , "oxygen", "mmHg", 1e5 , 0.01 );
    coarse_well.add_density( "glucose", "mmol", 1e5 , 10.0 );
    coarse_well.add_density( "lactate", "mmol", 1e5 , 0.0 );
    coarse_well.resize_space( 100, 1 , 1 );
    
    double dx = 100;
    coarse_well.resize_space_uniform( 0, 5000.0 , -dx/2.0 , dx/2.0 , -dx/2.0 , dx/2.0 , dx );
    std::vector<double> dirichlet_condition = { 1 , 1, 0 };
    
    int my_voxel_index = 0;
    coarse_well.add_dirichlet_node( my_voxel_index , dirichlet_condition );

    dirichlet_condition = { 0,0,0 };
    my_voxel_index = coarse_well.mesh.voxels.size()-1;
    coarse_well.add_dirichlet_node( my_voxel_index , dirichlet_condition );
    coarse_well.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;
    
    coarse_well.display_information( std::cout );
    
    
    for( int n = 0 ; n < coarse_well.mesh.voxels.size(); n++ )
    {
        std::cout << coarse_well(n)[0] << " ";
    }
    
    std::cout << std::endl;
    
    double t_temp = 0;
    while( t_temp < 100 )
      {
          coarse_well.simulate_diffusion_decay( 0.01 );
          t_temp += 0.01;
      }

      for( int n = 0 ; n < coarse_well.mesh.voxels.size(); n++ )

      {

            std::cout << coarse_well(n)[0] << " ";

      }
      std::cout << std::endl; 
	return; 
} 

void setup_tissue( void )
{
	// create some cells near the origin
	Cell* pCell;
    int retval;
    static Cell_Definition* pDefault = find_cell_definition("default");
    static Cell_Definition* pFibro = find_cell_definition("fibroblast");
    static Cell_Definition* pKRAS_pos = find_cell_definition("KRAS_positive");
    static Cell_Definition* pKRAS_neg = find_cell_definition("KRAS_negative");

    // std::cout << "setup_tissue:  *pFibro.type = " << *pFibro->type << std::endl;

    // 20,000 fibroblast seeding
    int kfib = 0;
    if (parameters.bools("fibroblast_seeding"))
    {
        std::cout << "creating fibroblasts" << std::endl;
        for (int i= -2666; i<2666; i+=16.82+20.8)  //rwh??
        {
            for (int j= -2666; j<2666; j+=16.82+20.8)  //rwh??
            {			
                kfib += 1;
                // pCell = create_cell(fibroblast);
                pCell = create_cell( *pFibro );
                pCell->assign_position(i,-500,j);
                pCell->functions.update_phenotype = tumor_energy_update_function_fibroblast;
            } 	
        }  
    }
    std::cout << "------ created # fibroblasts = " << kfib << std::endl;

    double cell_radius = cell_defaults.phenotype.geometry.radius; 
    // double initial_tumor_radius = 46; // parameters.doubles("initial_tumor_radius");
    double initial_tumor_radius = parameters.doubles("initial_tumor_radius");
    std::cout << "------ initial_tumor_radius = " << initial_tumor_radius << std::endl;

    //rwh
    // double number_of_organoid = 250; //parameters.doubles("number_of_organoid")
    // int number_of_organoid = 250; // parameters.doubles("number_of_organoid")
    int number_of_organoids = parameters.ints("number_of_organoids");
    std::cout << "------ number_of_organoids = " << number_of_organoids << std::endl;
    
    double xmin=1.e6;
    double ymin=1.e6;
    double zmin=1.e6;
    double xmax= -xmin;
    double ymax= -ymin;
    double zmax= -zmin;
	if (parameters.bools("organoid_cell_seeding"))
	{
            // std::cout << "creating CRCs" << std::endl;
            // rwh: if (parameters.doubles("organoid_cell_seeding_method") == 1)
            if (parameters.ints("organoid_cell_seeding_method") == 1)  // only KRAS_pos
            {
                std::cout << "  - seeding method = 1" << std::endl;
                //rwh for (int i = 0; i < number_of_organoid; i++) // seeding number of organoid cells specified in PhysiCell_settings.xml
                int kdx = 0;
                for (int idx = 0; idx < number_of_organoids; idx++) // seeding number of organoid cells specified in PhysiCell_settings.xml
			    {
                    //rwh - why??  and why in this loop?
                    std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
                    //std::cout << "creating " << positions.size() << " closely-packed organoid cells ... " << std::endl;
                    // create organoid
                        //rwh: create values in range: rand() % (max_number + 1 - minimum_number) + minimum_number; e.g.
                        // [2666, 5333+2666-1] = 2666, 
                        double xrand = (rand() % 5333) - 2666;
                        double yrand = (rand() % 961) - 480;
                        double zrand = (rand() % 5333) - 2666;
                        if (xrand < xmin) xmin = xrand;
                        if (xrand > xmax) xmax = xrand;
                        if (yrand < ymin) ymin = yrand;
                        if (yrand > ymax) ymax = yrand;
                        if (zrand < zmin) zmin = zrand;
                        if (zrand > zmax) zmax = zrand;
                    //std::cout << positions.size() << std::endl;
                    for( int i=0; i < positions.size(); i++ )
                    {
                        positions[i][0] += xrand;//(rand() % 5333) - 2666;
                        positions[i][1] += yrand;//(rand() % 961) - 480;
                        positions[i][2] += zrand;//(rand() % 5333) - 2666;
                        // pCell = create_cell(KRAS_positive);
                        pCell = create_cell( *pKRAS_pos );
                        pCell->assign_position( positions[i] );
                        pCell->functions.update_phenotype = tumor_energy_update_function_KRAS_positive;
                        kdx += 1;
                    }
			    }
                std::cout << " -----  # of KRAS_pos cells = " << kdx << std::endl;
            }
            //---------------------------------------------------
            // else if (parameters.doubles("organoid_cell_seeding_method") == 2)
            else if (parameters.ints("organoid_cell_seeding_method") == 2) // only KRAS_neg
            {
                std::cout << "  - seeding method = 2 (only KRAS_neg)" << std::endl;
                for (int i = 0; i < number_of_organoids; i++) // seeding number of organoid cells specified in PhysiCell_settings.xml
			    {
                
                    std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
                    //std::cout << "creating " << positions.size() << " closely-packed organoid cells ... " << std::endl;
                    // create organoid
                        double xrand = (rand() % 5333) - 2666;
                        double yrand = (rand() % 961) - 480;
                        double zrand = (rand() % 5333) - 2666;
                        if (xrand < xmin) xmin = xrand;
                        if (xrand > xmax) xmax = xrand;
                        if (yrand < ymin) ymin = yrand;
                        if (yrand > ymax) ymax = yrand;
                        if (zrand < zmin) zmin = zrand;
                        if (zrand > zmax) zmax = zrand;
                    //std::cout << positions.size() << std::endl;
                    for( int i=0; i < positions.size(); i++ )
                    {
                        positions[i][0] += xrand;//(rand() % 5333) - 2666;
                        positions[i][1] += yrand;//(rand() % 961) - 480;
                        positions[i][2] += zrand;//(rand() % 5333) - 2666;
                        // pCell = create_cell(KRAS_negative);
                        pCell = create_cell( *pKRAS_neg );
                        pCell->assign_position( positions[i] );
                        pCell->functions.update_phenotype = tumor_energy_update_function_KRAS_negative;
                    }
			    }
            }
            //---------------------------------------------------
            // else if (parameters.doubles("organoid_cell_seeding_method") == 3)
            else if (parameters.ints("organoid_cell_seeding_method") == 3) // both KRAS_pos, KRAS_neg
            {
                std::cout << "  - seeding method = 3 (both KRAS_pos, KRAS_neg)" << std::endl;

                int num_kras_pos_organoids = parameters.doubles("percent_KRAS_positive") * number_of_organoids;
                std::cout << "  num_kras_pos_organoids = " << num_kras_pos_organoids  << std::endl;

                for (int idx = 0; idx < number_of_organoids; idx++)
                {
                    double xrand = (rand() % 5333) - 2666;
                    double yrand = (rand() % 961) - 480;
                    double zrand = (rand() % 5333) - 2666;
                    std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
                    //std::cout << positions.size() << std::endl;
                    for( int i=0; i < positions.size(); i++ )
                    {
                        positions[i][0] += xrand;//(rand() % 5333) - 2666;
                        positions[i][1] += yrand;//(rand() % 961) - 480;
                        positions[i][2] += zrand;//(rand() % 5333) - 2666;

                        double num = UniformRandom();
                        
                        if (num <= parameters.doubles("percent_KRAS_positive"))
                        {
                            // pCell = create_cell(KRAS_negative);
                            pCell = create_cell( *pKRAS_pos );
                            pCell->assign_position( positions[i] );
                            pCell->functions.update_phenotype = tumor_energy_update_function_KRAS_positive;
                        }
                        else
                        {
                            pCell = create_cell( *pKRAS_neg );
                            pCell->assign_position( positions[i] );
                            pCell->functions.update_phenotype = tumor_energy_update_function_KRAS_negative;
                        }
                        
                        
                    }
                }
			}
        std::cout << "------- organoid seed cells ranges:\n";
        std::cout << " x: " << xmin << ", " << xmax << std::endl;
        std::cout << " y: " << ymin << ", " << ymax << std::endl;
        std::cout << " z: " << zmin << ", " << zmax << std::endl;
    }

        // pCell->custom_data[i_Lac_i] = pCell->phenotype.intracellular->get_double_parameter_value("Lactate");

        // phenotype.molecular.internalized_total_substrates[i_Lac] = pCell->custom_data[i_Lac_i]*cell_volume;
    //    freeRRCData (result);
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

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
}

void update_intracellular()
{
    // BioFVM Indices
    static int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" );
    static int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
	static int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" );
    static int lactate_substrate_index = microenvironment.find_density_index( "lactate");


    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        // Custom Data Indices
        static int i_Oxy_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_oxy" );
		static int i_Glc_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_glc" );
        static int i_Glu_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_glu" );
        static int i_Lac_i = (*all_cells)[i]->custom_data.find_variable_index( "intra_lac" );
        static int energy_vi = (*all_cells)[i]->custom_data.find_variable_index( "intra_energy" );

        if( (*all_cells)[i]->is_out_of_domain == false  )
        {
            // Cell Volume
            double cell_volume = (*all_cells)[i]->phenotype.volume.total;
            
            // Intracellular Concentrations
            double oxy_val_int = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]/cell_volume;
			double glc_val_int = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index]/cell_volume;
            double glu_val_int = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glutamine_substrate_index]/cell_volume;
            double lac_val_int = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index]/cell_volume;

            
            // Update SBML 
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Oxygen",oxy_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Glucose",glc_val_int);
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("Glutamine",glu_val_int);
            (*all_cells)[i]->phenotype.intracellular->set_parameter_value("Lactate",lac_val_int);
            
            
            // SBML Simulation
            (*all_cells)[i]->phenotype.intracellular->update();
            // Phenotype Simulation
            (*all_cells)[i]->phenotype.intracellular->update_phenotype_parameters((*all_cells)[i]->phenotype);
            
            
            // Internalized Chemical Update After SBML Simulation
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen") * cell_volume;
			(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glutamine_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glutamine") * cell_volume;
            (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate") * cell_volume;
            
            
            //Save custom data
            (*all_cells)[i]->custom_data[i_Oxy_i] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen");
            (*all_cells)[i]->custom_data[i_Glc_i] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glucose");
			(*all_cells)[i]->custom_data[i_Glu_i] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Glutamine");
            (*all_cells)[i]->custom_data[i_Lac_i] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Lactate");
            (*all_cells)[i]->custom_data[energy_vi] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy");

        }
    }
    
                
    //std::cout<<std::endl;
    //std::cout<<std::endl;
    //std::cout<<std::endl;
}
 

void tumor_energy_update_function_KRAS_positive( Cell* pCell, Phenotype& phenotype , double dt )
{
/* 	static int Start_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model);
	static int End_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model );

    double tr = phenotype.cycle.data.transition_rate( Start_index,End_index ); */
    //double i_Oxy_i = pCell->custom_data.find_variable_index( "oxygen_i_conc" );
    //std::cout << pCell->custom_data[i_Oxy_i] << std::endl;
    static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
/*     if( pCell->phenotype.death.dead == true )
    {
        pCell->functions.custom_cell_rule = NULL; 
		return; 
    } */
	
	
	// Apoptosis rate change is disabled.
  /*  static int lactate_index = microenvironment.find_density_index( "lactate" ); 
    double lactate_threshold = pCell->custom_data["lactate_threshold"];
    
    if( pCell->phenotype.death.dead == false && pCell->type == 3 )
    { 
        if( pCell->nearest_density_vector()[lactate_index] > lactate_threshold )
        {
            pCell->phenotype.death.rates[apoptosis_model_index] = 0.01;
        }
    } */
    
	return;
}

void tumor_energy_update_function_KRAS_negative( Cell* pCell, Phenotype& phenotype , double dt )
{
/* 	static int Start_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model);
	static int End_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model );

    double tr = phenotype.cycle.data.transition_rate( Start_index,End_index ); */
    //double i_Oxy_i = pCell->custom_data.find_variable_index( "oxygen_i_conc" );
    //std::cout << pCell->custom_data[i_Oxy_i] << std::endl;
    static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	
    // std::cout << "Kras Negative Apoptosis Rate = " << pCell->phenotype.death.rates[apoptosis_model_index] << std::endl;
    // std::cout << "Kras Negative Necrosis Rate = " << pCell->phenotype.death.rates[necrosis_model_index] << std::endl;
	
/*     if( pCell->phenotype.death.dead == true )
    {
        pCell->functions.custom_cell_rule = NULL; 
		return; 
    } */
	
	
	// Lactate-Apoptosis is disabled.
	
/* static int lactate_index = microenvironment.find_density_index( "lactate" ); 
    double lactate_threshold = pCell->custom_data["lactate_threshold"];
    
    if( pCell->phenotype.death.dead == false && pCell->type == 2 )
    { 
        if( pCell->nearest_density_vector()[lactate_index] > lactate_threshold )
        {
            // std::cout << "Dyiiinnggg KRAS_neg" << std::endl;
            pCell->phenotype.death.rates[apoptosis_model_index] = 0.01;
        }
    } */
    
	return;
}

void tumor_energy_update_function_fibroblast( Cell* pCell, Phenotype& phenotype , double dt )
{
/* 	static int Start_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model);
	static int End_index = live.find_phase_index( PhysiCell_constants::live_cells_cycle_model );

    double tr = phenotype.cycle.data.transition_rate( Start_index,End_index ); */
    //double i_Oxy_i = pCell->custom_data.find_variable_index( "oxygen_i_conc" );
    //std::cout << pCell->custom_data[i_Oxy_i] << std::endl;
    static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
/*     if( pCell->phenotype.death.dead == true )
    {
        pCell->functions.custom_cell_rule = NULL; 
		return; 
    } */
	
	
    // Apoptosis rate change is disabled.
	/*
    static int lactate_index = microenvironment.find_density_index( "lactate" ); 
    double lactate_threshold = pCell->custom_data["lactate_threshold"];
    
    if( pCell->phenotype.death.dead == false && pCell->type == 1 )
    { 
        if( pCell->nearest_density_vector()[lactate_index] > lactate_threshold )
        {
            // std::cout << "Dyiiinnggg fibro" << std::endl;
            pCell->phenotype.death.rates[apoptosis_model_index] = 0.01;
        }
    } */
    
	return;
}
