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

std::vector<float3> TriangularLattice(float terrain_rad, float world_size, int baseCount, int numLayers) {
    std::vector<float3> positions;
    float diameter = 2 * terrain_rad;

    //1% safety margin before generation
    float safetyMargin = 1.01;
    float inLayerSpacing = diameter * safetyMargin;
    float verticalSpacing = sqrt(3)/2 * diameter * safetyMargin;
    float layerHeight = 0; 
    bool evenLayer = true;

    for (int layer = 0; layer < numLayers; ++layer) {
        int particlesInLayer = evenLayer ? baseCount : (baseCount - 1);
        for (int i = 0; i <particlesInLayer; ++i) {
            //calculating the x position
            float xOffset = evenLayer ? (i - baseCount / 2.0f) * inLayerSpacing : ((i + 0.5) - baseCount / 2.0f) * inLayerSpacing;
            positions.push_back(make_float3(xOffset, 0, layerHeight));
        }
        layerHeight += verticalSpacing;
        evenLayer = !evenLayer;
    }

    return positions;
}

int main() {
    DEMSolver DEMSim;
    // Output less info at initialization
    DEMSim.SetVerbosity("INFO");
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
                                      {0, 30.0*world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Generate terrain particles
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.0e3 * 4 / 3 * PI,
                                                      terrain_rad, mat_type_terrain); 

    auto latticePositions = TriangularLattice(terrain_rad, world_size, 61, 15)
   
    //add clumps                                               
    DEMSim.AddClumps(template_terrain, latticePositions);

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
