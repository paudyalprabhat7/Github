//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

///////////////////////////////////////////////////////////////////////////////////////////////////////
void write_csv(std::string filename, std::vector<std::pair<std::string, std::vector<float>>> dataset){
    // Make a CSV file with one or more columns of integer values
    // Each column of data is represented by the pair <column name, column data>
    //   as std::pair<std::string, std::vector<int>>
    // The dataset is represented as a vector of these columns
    // Note that all columns should be the same size
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < dataset.size(); ++j)
    {
        myFile << dataset.at(j).first;
        if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << "\n";
    
    // Send data to the stream
    for(int i = 0; i < dataset.at(0).second.size(); ++i)
    {
        for(int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).second.at(i);
            if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }
    
    // Close the file
    myFile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// | OUTPUT_CONTENT:: VEL //|  OUTPUT_CONTENT:: ACC 
// | NORMAL


int main() {

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetErrorOutVelocity(10000.);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV  |  OUTPUT_CONTENT:: ABS_ACC );
    DEMSim.SetContactOutputContent(OWNER | FORCE | POINT | COMPONENT  );
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.EnsureKernelErrMsgLineNum();

    // E, nu, CoR, mu, Crr...
    auto mat_type_screw = DEMSim.LoadMaterial({{"E", 4.1e9}, {"nu", 0.394}, {"CoR", 0.3}, {"mu", 0.04}, {"Crr", 0.01}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 100e9}, {"nu", 0.45}, {"CoR", 0.5}, {"mu", 0.4}, {"Crr", 0.175}});
    auto mat_type_container = DEMSim.LoadMaterial({{"E", 100e9}, {"nu", 0.45}, {"CoR", 0.5}, {"mu", 0.4}, {"Crr", 0.175}});
    // If you don't have this line, then CoR,Crr,mu between screw material and granular material will be average of the two.
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_screw, mat_type_terrain, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_screw, mat_type_terrain, 0.01);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_screw, mat_type_terrain, 0.04);

    DEMSim.SetMaterialPropertyPair("CoR", mat_type_container, mat_type_terrain, 0.36);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_container, mat_type_terrain, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_container, mat_type_terrain, 1.00);
    
    // 1- Boundary:
    float step_size = 1e-6;
    //float Sim_CED = 0; //Pa

    //lateral and vertical extents
    float terrain_rad = 0.005;

    double xBoundary = 122.0 * terrain_rad;
    double yBoundary = 12.0 * terrain_rad;
    double zBoundary = 30.0 * terrain_rad;
    double safety_margin= 0.05;

    DEMSim.InstructBoxDomainDimension(xBoundary, yBoundary, zBoundary + safety_margin);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_container);
    auto bot_wall = DEMSim.AddBCPlane(make_float3(0,0,-zBoundary/2.0),make_float3(0,0,1),mat_type_container);

    //2-Bottom Wall:
    bot_wall->SetMass(1000);//1000 is the assumed mass of boundary. 
    bot_wall->SetMOI(make_float3(1./12.*1000.0*yBoundary*yBoundary, 1./12.*1000.0*xBoundary*xBoundary, 1./12.*1000.0*(xBoundary*xBoundary+yBoundary*yBoundary)));
    auto bot_wall_tracker = DEMSim.Track(bot_wall);

    //load terrain
    
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.6e3 * 4 / 3 * 3.14,
                                                  terrain_rad, mat_type_terrain);

    // Define sampling volume dimensions within the container
    float sample_halfheight = zBoundary / 4; // Example: Fill the bottom quarter of the container height with particles
    float3 sample_center = make_float3(0, 0, -zBoundary / 2 + sample_halfheight + 0.05); // Adjust center based on bottom wall position and desired fill height
    float sample_halfwidth_x = xBoundary / 2 * 0.95; // Slightly smaller than container to avoid wall overlap
    float sample_halfwidth_y = yBoundary / 2 * 0.95; // Slightly smaller than container to avoid wall overlap

    // Use DEMBoxHCPSampler to sample within the defined volume
    auto input_xyz = DEMBoxHCPSampler(sample_center, make_float3(sample_halfwidth_x, sample_halfwidth_y, sample_halfheight), 2.01 * terrain_rad);

    // Add sampled clumps to the simulation using your defined template
    DEMSim.AddClumps(template_terrain, input_xyz);
    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;
    

    // A custom force model can be read in through a file and used by the simulation. Magic, right?
    //auto my_force_model = DEMSim.ReadContactForceModel("SJKR.cu");
    // This custom force model still uses contact history arrays, so let's define it
    //my_force_model->SetPerContactWildcards({"delta_tan_x", "delta_tan_y", "delta_tan_z", "delta_time"});
    //std::cout << "Force Module loaded:" << n_particles << std::endl;
    // 4- Compressor:
    // Now add a plane to compress the sample
    //auto compressor = DEMSim.AddExternalObject();
    //compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, 1), mat_type_container);
    //compressor->SetFamily(10);
    //DEMSim.DisableContactBetweenFamilies(0, 10);
    //DEMSim.DisableContactBetweenFamilies(1, 10);
    //DEMSim.DisableContactBetweenFamilies(2, 10);
    //DEMSim.DisableContactBetweenFamilies(3, 10);
    //DEMSim.SetFamilyFixed(10);
    //auto compressor_tracker = DEMSim.Track(compressor);


    // 5- Simulation initilization:
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    //DEMSim.SetMaxVelocity(15.);
    DEMSim.SetExpandSafetyAdder(5.);

    DEMSim.Initialize();

    // 6- Drop Particles:
    std::cout << "//////////////// Particles generating //////////////////: " << std::endl; 

    path out_dir = current_path();
    out_dir += "/Containerwithbottomplane";
    create_directory(out_dir);

    float settle_time = 2.0;
    unsigned int fps = 20;
    float frame_time = 1.0/fps;

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;

    for(float t=0; t<settle_time; t+=frame_time){
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200];
        sprintf(filename, "%s/DEM_particlelattice_out_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEM_particlelattice_%04d.vtk", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);

    }
    
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();


    DEMSim.ShowTimingStats();
  
    return 0;

}
