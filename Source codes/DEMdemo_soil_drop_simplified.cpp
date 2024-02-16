//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This demo tries to show the strain distribution in the granular material when
// affected by a compressor.
// =============================================================================

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


int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);

    path out_dir = current_path();
    out_dir += "/DemoOutput_Indentation_boxdrop";
    create_directory(out_dir);

    // E, nu, CoR, mu, Crr...
    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.1}, {"mu", 0.4}, {"Crr", 0.0}});
    auto mat_type_granular_1 = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.1}, {"mu", 0.45}, {"Crr", 0.5}});
    auto mat_type_granular_2 = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.1}, {"mu", 0.45}, {"Crr", 0.5}});
    // CoR is a pair-wise property, so it should be mentioned here
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_cube, mat_type_granular_1, 0.8);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_cube, mat_type_granular_2, 0.8);

    float granular_rad = 0.005;  // 0.002;
    auto template_granular = DEMSim.LoadSphereType(granular_rad * granular_rad * granular_rad * 2.6e3 * 4 / 3 * 3.14,
                                                   granular_rad, mat_type_granular_1);

    float step_size = 1e-5;
    const double world_size = 0.15;
    const float fill_height = 0.075;
    const float chamber_bottom = -world_size / 2.;
    const float fill_bottom = chamber_bottom + granular_rad;

    DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_granular_2);

    // define the upper cube
    auto cube_upper = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cube.obj"), mat_type_cube);
    std::cout << "Total num of triangles in upper cube: " << cube_upper->GetNumTriangles() << std::endl;
    // Make the cube about 2cm by 1cm
    float cube_width = 0.02;
    float cube_height = 0.01;
    //double cube_speed = 0.25;  // 0.1 and 0.02, try them too... very similar though
    cube_upper->Scale(make_float3(cube_width, cube_width, cube_height));
    cube_upper->SetFamily(3);
    DEMSim.SetFamilyFixed(3);

    //DEMSim.SetFamilyPrescribedLinVel(11, "0", "0", to_string_with_precision(-cube_speed));
    // Track the cube
    auto cube_upper_tracker = DEMSim.Track(cube_upper);

    //define the lower cube
    auto cube_lower = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cube.obj"), mat_type_cube);
    std::cout << "Total num of triangles in lower cube: " << cube_lower->GetNumTriangles() << std::endl;
    float cube_width_l = 0.02;
    float cube_height_l = 0.01;
    //double cube_speed = 0.25;  // 0.1 and 0.02, try them too... very similar though
    cube_lower->Scale(make_float3(cube_width, cube_width, cube_height));
    cube_lower->SetFamily(4);
    DEMSim.SetFamilyFixed(4);

    //track the lower cube
    auto cube_lower_tracker = DEMSim.Track(cube_lower); 

    HCPSampler sampler(3.f * granular_rad);
    float3 fill_center = make_float3(0, 0, fill_bottom + fill_height / 2);
    const float fill_radius = world_size / 2 - 2. * granular_rad;
    auto input_xyz = sampler.SampleCylinderZ(fill_center, fill_radius, fill_height);

    auto particles = DEMSim.AddClumps(template_granular, input_xyz);
    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;

    // After creating particles
    // Declare and initialize the particle tracker and related variables
    auto particle_tracker = DEMSim.Track(particles);
    unsigned int num_particles = input_xyz.size();
    particles->SetFamily(1);

    // Position the upper cube
    float upper_cube_z_position = fill_bottom + fill_height + 0.01 + cube_height / 2.0;
    cube_upper_tracker->SetPos(make_float3(0, 0, upper_cube_z_position));

    // Position the lower cube
    float depth_80_percent = fill_bottom + 0.8 * fill_height - cube_height_l / 2.0;
    cube_lower_tracker->SetPos(make_float3(0, 0, depth_80_percent));

    DEMSim.DisableContactBetweenFamilies(1,3);

    //gravity and others
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetCDUpdateFreq(20);
    // You usually don't have to worry about initial bin size. But sometimes if you can set the init bin size so that
    // the kT--dT work at a sweet collaboration pattern, it could make the solver run faster.
    DEMSim.SetInitBinNumTarget(5e7);
    DEMSim.Initialize();


    //letting the first cube drop
    DEMSim.ChangeFamily(3,5);
    
    const double total_simulation_time = 5.0; // Total time for the simulation
    const unsigned int fps = 200;             // Frames per second for output
    const float frame_time = 1.0 / fps;       // Time interval for outputs
    unsigned int out_steps = static_cast<unsigned int>(1.0 / (fps * step_size));

    std::cout << "Starting the simulation..." << std::endl;

    // Simplified Simulation loop
    for (double t = 0; t < total_simulation_time; t += step_size) {
        DEMSim.DoDynamics(step_size);

        // Output data at specified intervals
        if (static_cast<unsigned int>(t / step_size) % out_steps == 0) {
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << static_cast<unsigned int>(t / frame_time) << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), static_cast<unsigned int>(t / frame_time));
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), static_cast<unsigned int>(t / frame_time));
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));
        }
    }

    std::cout << "Simulation completed." << std::endl;
    return 0;
}
