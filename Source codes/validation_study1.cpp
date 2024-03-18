//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// A demo that is basically the hello-world script for DEME.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>
#include <cmath>

#include <filesystem>
#include <cstdio>
#include <time.h>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"OWNER", "FORCE", "POINT", "COMPONENT", "NORMAL", "TORQUE"});
    DEMSim.EnsureKernelErrMsgLineNum();
    DEMSim.SetErrorOutVelocity(2000.);

    // srand(time(NULL));
    srand(4150);

    //defining the material type
    auto mat_type_1 =
        DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.3}, {"Crr", 0.01}});
    
    
    //definition of the sphere type
    auto sph_type_1 = DEMSim.LoadSphereType(11728., 1., mat_type_1);

    // Adjusted parameters for custom lattice stacking
    const float R = 0.05; // Particle radius
    const float particleDiameter = 2 * R;
    const int baseLayerCountX = 60; // specifying the number of particles in the bottom layer x-direction
    const int totalLayers = 15; // Total number of layers to stack

    std::vector<float3> positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Generate the custom lattice with additional spacing to avoid incorrect contact force resolution
    for (int layer = 0; layer < totalLayers; ++layer) {
        // Alternate between baseLayerCountX and baseLayerCountX - 1 based on whether the layer index is even or odd
        int layerCountX = baseLayerCountX - (layer % 2);
        // Calculation for layerCountY remains the same
        int layerCountY = static_cast<int>(std::round(0.1 * layerCountX));

        // Calculate an offset for odd layers to create a densely packed structure
        float layerOffset = 0.0f;
        if (layer % 2 == 1) {
            // Offset odd layers by half the particle diameter to create a hexagonal close packing
            layerOffset = 0.5 * particleDiameter;
        }

        // Additional space of 10% of the particle diameter for separation
        float additionalSpace = 0.1 * particleDiameter;

        for (int y = 0; y < layerCountY; ++y) {
            for (int x = 0; x < layerCountX; ++x) {
                // Adjust xPos calculation by adding the layerOffset and consider the additional spacing
                float xPos = (x * particleDiameter + layerOffset) + additionalSpace;
                // Incorporate the additional spacing in yPos and zPos calculations
                float yPos = (y * sqrt(3) * R) + additionalSpace;
                // Adjust zPos for additional spacing between layers
                float zPos = (layer * particleDiameter) + (layer > 0 ? additionalSpace : 0.0);

                positions.push_back(make_float3(xPos, yPos, zPos));
                clump_types.push_back(sph_type_1);
            }
        }
    }

    auto particles = DEMSim.AddClumps(clump_types, positions);
    
    DEMSim.SetInitTimeStep(2e-5);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.8));
    DEMSim.SetCDUpdateFreq(10);
    DEMSim.SetMaxVelocity(6.);
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.SetExpandSafetyMultiplier(1.2);
    DEMSim.SetIntegrator("centered_difference");

    DEMSim.Initialize();
    
    path out_dir = current_path();
    out_dir += "/DEM_particlelattice";
    create_directory(out_dir);

    float settle_time = 2.0;
    unsigned int fps = 20;
    float frame_time = 1.0 / fps;

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;

    //let's settle the lattice
    for(float t=0; t < settle_time; t += frame_time) {
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200];
        sprintf(filename, "%s/DEM_particlelattice_out_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEM_particlelattice_%04d.vtk", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }
    
}