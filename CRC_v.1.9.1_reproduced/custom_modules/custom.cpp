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
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;


#include "rrc_api.h"
#include "rrc_types.h"
// #include "rrc_utilities.h"
extern "C" rrc::RRHandle createRRInstance();

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	build_cell_definitions_maps(); 
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

   // create some cells near the origin
	Cell* pCell;
    static Cell_Definition* pDefault = find_cell_definition("default");
    static Cell_Definition* pFibro = find_cell_definition("fibroblast");
    static Cell_Definition* pKRAS_pos = find_cell_definition("KRAS_positive");
    static Cell_Definition* pKRAS_neg = find_cell_definition("KRAS_negative");


    // Fibroblast seeding
    int kfib = 0;
    if (parameters.bools("fibroblast_seeding"))
    {
        std::cout << "creating fibroblasts" << std::endl;
        for (int i= -2666; i<2666; i+=16.82+20.8)
        {
            for (int j= -2666; j<2666; j+=16.82+20.8)
            {			
                kfib += 1;
                pCell = create_cell( *pFibro );
                pCell->assign_position(i,-500,j);
                if (parameters.bools("add_SBML_Modeling") == true)
                {
                    pCell->phenotype.intracellular->start(); // INTRACELLULAR
                }
                pCell->functions.update_phenotype = tumor_energy_update_function_fibroblast;
            } 	
        }  
    }
    std::cout << "------ created # fibroblasts = " << kfib << std::endl;

    double cell_radius = cell_defaults.phenotype.geometry.radius; 
    double initial_tumor_radius = parameters.doubles("initial_tumor_radius");
    std::cout << "------ initial_tumor_radius = " << initial_tumor_radius << std::endl;


    int number_of_organoids = parameters.ints("number_of_organoids");
    std::cout << "------ number_of_organoids = " << number_of_organoids << std::endl;
    
    double xmin=1.e6; // Work here @Furkan 
    double ymin=1.e6;
    double zmin=1.e6;
    double xmax= -xmin;
    double ymax= -ymin;
    double zmax= -zmin;
	if (parameters.bools("organoid_cell_seeding"))
	{
            // std::cout << "creating CRCs" << std::endl;
            if (parameters.ints("organoid_cell_seeding_method") == 1)  // only KRAS_pos
            {
                std::cout << "  - seeding method = 1" << std::endl;

                int kdx = 0;
                for (int idx = 0; idx < number_of_organoids; idx++) // seeding number of organoid cells specified in PhysiCell_settings.xml
			    {
                    std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 

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
                        if (parameters.bools("add_SBML_Modeling") == true)
                        {
                            pCell->phenotype.intracellular->start(); // INTRACELLULAR
                        }
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
                        if (parameters.bools("add_SBML_Modeling") == true)
                        {
                            pCell->phenotype.intracellular->start(); // INTRACELLULAR
                        }
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
                            if (parameters.bools("add_SBML_Modeling") == true)
                            {
                                pCell->phenotype.intracellular->start(); // INTRACELLULAR
                            }
                            pCell->functions.update_phenotype = tumor_energy_update_function_KRAS_positive;
                        }
                        else
                        {
                            pCell = create_cell( *pKRAS_neg );
                            pCell->assign_position( positions[i] );
                            if (parameters.bools("add_SBML_Modeling") == true)
                            {
                                pCell->phenotype.intracellular->start(); // INTRACELLULAR
                            }
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

    return;
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
    return;
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
