//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Enter code description here
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstudio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

int main() {
    DEMSolver DEMSim; //declare and initialize the object DEMSim of the DEMSolver class
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV); 
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    //DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.EnsureKernelErrMsgLineNum(); 

    //load the material properties
    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 2.1e5}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, "Crr", 0.01});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e5}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.65}, "Crr", 0.01});
    auto mat_type_analyticalb = DEMSim.LoadMaterial({{"E", 2.1e5}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, "Crr", 0.01});

    //step size
    float step_size = 1e-4;
    float world_size = 2;

    //Analytical boundary definition
    DEMSim.InstructBoxDomainDimension({0, world_size}, {0, world_size}, {0, world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_analyticalb);

    //terrain definition
    float terrain_rad = 0.08;
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.69e3 * 4/3 * 3.141, terrain_rad, mat_type_terrain);

    //terrain sampling
    float sample_halfheight = world_size/4;
    float3 sample_center = (world_size/2, world_size/2, sample_halfheight);
    float sample_halfwidth = world_size/2;
    auto input_xyz = DEMBoxHCPSampler(sample_center, make_float3(sample_halfwidth, sample_halfwidth, sample_halfheight), 2.01*terrain_rad);
    DEMSim.AddClumps(template_terrain, input_xyz);

    //initialization of simulation
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(0, 0 , -9.81);
    DEMSim.SetMaxVelocity(15.);

    DEMSim.Initialize();

    //creating the output directory
    path out_dir = current_path();
    out_dir += "/Output_settlingphase";
    create_directory(out_dir);

    //visualization frame time
    float settle_time = 2.0;
    unsigned int fps = 24;
    float frame_time = 1.0 / fps;
    std::cout << "Output at " << fps <<" FPS" << std::endl;
    
    unsigned int curr_frame = 0;

    //main loop for settling
    for (float t = 0; t<settle_time; t+=frame_time) {
        char filename[200];
        sprintf(filename, "%s/settlingphase_output_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        curr_frame++;

        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }

    //post settling housekeeping
    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();
    std::cout << "DEMdemo_BallDrop exiting..." << std::endl;
    return 0;

}