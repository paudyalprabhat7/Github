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

std::vector<float3> Generate2DTriangularLattice(float particleRadius, int layers) {
    std::vector<float3> positions;
    float diameter = 2.0f * particleRadius;
    // Apply a 5% safety margin
    float spacing = diameter * 1.05f; //changing the spacing 
    float layerHeight = spacing * sqrt(3) / 2; // Vertical spacing for a triangular lattice

    for (int layer = 0; layer < layers; ++layer) {
        int particlesInLayer = layer % 2 == 0 ? 61 : 60;
        // Determine the offset for this layer to stagger the particles
        float offsetX = (layer % 2 == 0) ? 0 : spacing / 2;

        for (int particle = 0; particle < particlesInLayer; ++particle) {
            // Calculate positions, keeping Y constant as it's a 2D lattice
            float x = particle * spacing + offsetX;
            float z = layer * layerHeight; // Stack layers vertically
            float3 position;
            position.x = x;
            position.y = 0.0f;
            position.z = z;
            positions.push_back(position);
        }
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
    out_dir += "/DemoOutput_ParticleSettle_customlattice";
    create_directory(out_dir);

    // Material properties for the terrain
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});

    float terrain_rad = 0.003;
    int layers = 15;

    float diameter = 2.0f * terrain_rad;
    float spacing = diameter * 1.05f;
    float layerHeight = spacing * sqrt(3) / 2;

    //calculating the extents
    int particlesInLayer = 61;
    float maxXExtent = (particlesInLayer - 1) * spacing;

    //calculating the z extent
    float maxZExtent = (layers -1) *layerHeight;

    float step_size = 2e-6;
    
    float buffer = terrain_rad;
    DEMSim.InstructBoxDomainDimension({-buffer, maxXExtent + buffer}, {-buffer, maxXExtent + buffer}, {0, maxZExtent + buffer});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Generate terrain particles
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.0e3 * 4 / 3 * PI,
                                                      terrain_rad, mat_type_terrain); 

    auto latticePositions = Generate2DTriangularLattice(terrain_rad, layers);
   
    //add clumps                                               
    DEMSim.AddClumps(template_terrain, latticePositions);

    //output initial positions for verification
    for (const auto& pos : latticePositions) {
        std::cout << "Particle Position: X = " << pos.x << ", Y = " << pos.y << ", Z = " << pos.z << std::endl;
    }

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

    std::cout << "Total num of particles: " << latticePositions.size() << std::endl;
    std::cout << "Particle settling simulation completed." << std::endl;
    return 0;
}
