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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom(); 
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment();
	
    // START -------- 1D Microenvironment ------------- START //
    
    Microenvironment coarse_well;
    Microenvironment old_coarse_well;
    Microenvironment new_coarse_well;
    
    coarse_well.name = "coarse well";
    coarse_well.spatial_units = "micron";
    coarse_well.mesh.units = "micron";
    coarse_well.time_units = "min";

    coarse_well.set_density( 0 , "oxygen", "mmHg", 1e5 , 0.00 );
    coarse_well.add_density( "glucose", "mM", 3e4 , 0.0 );
    coarse_well.add_density( "glutamine", "mM", 3e4 , 0.0);
    coarse_well.add_density( "lactate", "mM", 3e4 , 0.0 );
    coarse_well.resize_space( 100, 1 , 1 );
    
    double dx = 32;
    coarse_well.resize_space_uniform( -512.0, 10208.0 , -dx/2.0 , dx/2.0 , -dx/2.0 , dx/2.0 , dx );
    std::vector<double> dirichlet_condition = { 38 , 0, 0, 0 };
    
    int my_voxel_index = 319;
    coarse_well.set_substrate_dirichlet_activation(1,false);
    coarse_well.set_substrate_dirichlet_activation(2,false);
    coarse_well.set_substrate_dirichlet_activation(3,false);

    
    my_voxel_index = coarse_well.mesh.voxels.size()-1;
    coarse_well.add_dirichlet_node( my_voxel_index , dirichlet_condition );
    coarse_well.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;
    
    for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
    {
        coarse_well(m)[0]=38; // oxygen
        coarse_well(m)[1]=16.897255; // glucose
        coarse_well(m)[2]=5.40;
        coarse_well(m)[3]=0; // lactate
        //std::cout<< "turned" << std::endl;
    }
    
    coarse_well.display_information( std::cout );
    coarse_well.write_to_matlab("output/output00000000_microenvironment1.mat");
    
    // END -------- 1D Microenvironment ------------- END //
    
    //First sourcing (This is a special projection in the initial time step)
    for( int n = 0; n < coarse_well.mesh.voxels.size() ; n++ )
    {
       microenvironment(n)[0]=coarse_well(n)[0]; // oxygen
       microenvironment(n)[1]=coarse_well(n)[1]; // glucose
       microenvironment(n)[2]=coarse_well(n)[2]; // glutamine
       microenvironment(n)[3]=coarse_well(n)[3]; // lactate
    }    
    
    //setup_1D_microenvironment();
    //update_coarse_microenvironment();
    
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 32; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true );   
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	display_citations(); 
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
            
/*             for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
            {
                double mic_cen = microenvironment.mesh.voxels[n].center[1];
                for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
                {
                    double cmic_cen = coarse_well.mesh.voxels[m].center[0];
                    if (cmic_cen == mic_cen)
                    {
                       microenvironment(n)[0]=coarse_well(m)[0]; // oxygen
                       microenvironment(n)[1]=coarse_well(m)[1]; // glucose
                       microenvironment(n)[2]=coarse_well(m)[2]; // lactate
                       //counter+=15;
                    }
                }
            } */
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
  
                // START ------  CMicEnv Saving --------- START //
                //std::cout << "oxygen concentrations:" << std::endl; 
                /* for( int n = 0 ; n < coarse_well.mesh.voxels.size(); n++ )
                {
                std::cout << coarse_well(n)[0] << " ";
                }
                std::cout << std::endl; */

                sprintf( filename , "%s/output%08u_microenvironment1.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
                
                coarse_well.write_to_matlab(filename);
                // END ------  CMicEnv Saving --------- END //

			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}
            
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
            
            //Coarsening          
            if (fabs( PhysiCell_globals.current_time - 4  ) < 0.0000000001)
            {
            for( int n = -496; n < 497 ; n = n + 32 )////////
            {
                std::vector<double> v = {0, 0, 0, 0};
                for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
                {
                    double mic_cen_y = microenvironment.mesh.voxels[m].center[1];
                    if (mic_cen_y == n)
                    {
                        v[0]+=microenvironment(m)[0]*microenvironment.mesh.voxels[m].volume;
                        v[1]+=microenvironment(m)[1]*microenvironment.mesh.voxels[m].volume;
                        v[2]+=microenvironment(m)[2]*microenvironment.mesh.voxels[m].volume;
                        v[3]+=microenvironment(m)[3]*microenvironment.mesh.voxels[m].volume;
                    }
                }
                //std::cout << v[0] << std::endl;
                v[0]=v[0]/coarse_well.mesh.voxels[1].volume;
                v[1]=v[2]/coarse_well.mesh.voxels[1].volume;
                v[2]=v[2]/coarse_well.mesh.voxels[1].volume;
                v[2]=v[2]/coarse_well.mesh.voxels[1].volume;
                //std::cout << v[0] << std::endl;
                
                for ( int m = 0; m < 32 ; m++)
                {
                    coarse_well(m)[0]=v[0]; // oxygen
                    coarse_well(m)[1]=v[1]; // glucose
                    coarse_well(m)[2]=v[2]; // glutamine
                    coarse_well(m)[3]=v[3]; // lactate
                }
            }

            }
            

            coarse_well.simulate_diffusion_decay( diffusion_dt ); 



            // Projection---
           

            
            //std::cout<< PhysiCell_globals.current_time << std::endl;;
            
            //int counter = 15;
            // update the microenvironment according to coarse microenvironment
/*             if (fabs( PhysiCell_globals.current_time - 15  ) < 0.0000000001)
            {
                for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
                {
                    double mic_cen = microenvironment.mesh.voxels[n].center[1];
                    for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
                    {
                        double cmic_cen = coarse_well.mesh.voxels[m].center[0];
                        if (cmic_cen == mic_cen)
                        {
                           microenvironment(n)[0]=coarse_well(m)[0]; // oxygen
                           microenvironment(n)[1]=coarse_well(m)[1]; // glucose
                           microenvironment(n)[2]=coarse_well(m)[2]; // lactate
                           //counter+=15;
                        }
                    }
                }
            } */
            
             //
/*              // Start-- Refreshment
             if (fabs( PhysiCell_globals.current_time - 1728  ) < 0.0000000001)
            {
            for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
            {
                coarse_well(m)[0]=38; // oxygen
                coarse_well(m)[1]=16.897255; // glucose
                coarse_well(m)[2]=5.40; //glutamine
                coarse_well(m)[3]=0; // lactate
            }
            } */

             // End-- Refreshment
             
            
			// std::cout << "where is the bug?" << std::endl;
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
            //simulate_SBML_for_all_cells();
			
            
			/*
			  Custom add-ons could potentially go here. 
			*/			
			
			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}
