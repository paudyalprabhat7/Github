//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// ==========================================================================================
// 2d lattice settle before compression phase and load application using in-built box sampler 
// ==========================================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <random>

using namespace deme;
using namespace std::filesystem;

int main() {
    DEMSolver DEMSim;
    // Output less info at initialization
    DEMSim.SetVerbosity("ERROR");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetMeshOutputFormat("VTK");
    DEMSim.SetErrorOutVelocity(20000.);

    path out_dir = current_path();
    out_dir += "/DemoOutput_ParticleSettle";
    create_directory(out_dir);

    // Material properties for the terrain
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});

    double terrain_rad = 0.006 / 2.;

    float step_size = 2e-7;
    double world_size = terrain_rad * 124.0;
    DEMSim.InstructBoxDomainDimension({-world_size / 2., world_size / 2.}, {-world_size / 2., world_size / 2.},
                                      {0, 30.0*terrain_rad});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Generate terrain particles
    auto templates_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.0e3 * 4 / 3 * PI,
                                                      terrain_rad, mat_type_terrain); 

    HCPSampler sampler(2.01 * terrain_rad);
    float sample_halfwidth = world_size / 2 - 2 * terrain_rad;
    float fullheight = world_size * 6.;
    auto sample_center = make_float3(0, 0, fullheight / 2 + 1 * terrain_rad);
    auto input_xyz = sampler.SampleBox(sample_center, make_float3(sample_halfwidth, 0.f, fullheight / 2.));

    // Uniform selection of templates for each particle
    DEMSim.LoadSphereType(uniform_terrain_rad * uniform_terrain_rad * uniform_terrain_rad * 2.0e3 * 4 / 3 * PI,
                                                      uniform_terrain_rad, mat_type_terrain);

    std::vector<std::shared_ptr<DEMClumpTemplate>> template_to_use(input_xyz.size(), terrain_rad);

    //add clumps                                               
    DEMSim.AddClumps(template_to_use, input_xyz);

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.Initialize();

    float settle_time = 1.0;
    unsigned int fps = 20;
    float frame_time = 1.0 / fps;
    unsigned int out_steps = static_cast<unsigned int>(1.0 / (fps * step_size));
    unsigned int currframe = 0;

    // Simulation loop for settling
    for (float t = 0; t < settle_time; t += frame_time) {
        char filename[200];
        sprintf(filename, "%s/ParticleSettle_output_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        currframe++;
        DEMSim.DoDynamicsThenSync(frame_time);
    }

    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;
    std::cout << "Particle settling simulation completed." << std::endl;
    return 0;
}
