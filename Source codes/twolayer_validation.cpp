//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Lattice generation custom
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

    // srand(time(NULL));
    srand(4150);

    // Special material: has a cohesion param
    auto mat_type_1 =
        DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.3}, {"Crr", 0.01}});
    auto mat_type_2 =
        DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.4}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
    std::shared_ptr<DEMMaterial> mat_type_3 = DEMSim.Duplicate(mat_type_2);
    // If you don't have this line, then CoR between thw 2 materials will take average when they are in contact
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_1, mat_type_2, 0.6);
    // Even though set elsewhere to be 50, this pairwise value takes precedence when two different materials are in
    // contact.
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_1, mat_type_2, 0.5);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_1, mat_type_3, 0.5);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_2, mat_type_3, 0.5);

    
    //definition of the sphere type
    auto sph_type_1 = DEMSim.LoadSphereType(11728., 1., mat_type_1);

    //simulation misc
    DEMSim.SetInitTimeStep(2e-5);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.8));
    DEMSim.SetCDUpdateFreq(10);
    DEMSim.SetMaxVelocity(6.);
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.SetExpandSafetyMultiplier(1.2);
    DEMSim.SetIntegrator("centered_difference");
    
    //settlement parameters to add to the loop
    unsigned int globalFrame = 0;
    float settle_time_per_layer = 2.0;
    unsigned int fps = 20;
    float frame_time = 1.0/fps;

    // Adjusted parameters for custom lattice stacking
    const float R = 0.05; // Particle radius
    const float particleDiameter = 2 * R;
    const int baseLayerCountX = 16; // Starting with 16 particles in the bottom layer
    const int totalLayers = 2;

    std::vector<float3> positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    path out_dir = current_path() / "DEM_particlelattice_simulation";
    create_directory(out_dir); // Ensure this directory exists

    DEMSim.Initialize();
    
    for (int layer = 0; layer < totalLayers; ++layer) {
        positions.clear();
        clump_types.clear();

        int layerCountX = baseLayerCountX - layer;

        for (int x = 0; x < layerCountX; ++x) {
            float xPos = x * particleDiameter;
            float yPos = 0; // Assuming a single row for simplicity, adjust as needed
            float zPos = layer * particleDiameter; // Stack the layers directly on top of each other

            positions.push_back(make_float3(xPos, yPos, zPos));
            clump_types.push_back(sph_type_1); // Assuming sph_type_1 is defined as before
        }

        //add layer n to the simulation
        auto particles = DEMSim.AddClumps(clump_types, positions);

        // Simulate and output for the current layer
        for (float t = 0; t < settle_time_per_layer; t += frame_time) {
            char filename[200], meshfilename[200];
            sprintf(filename, "%s/DEM_particlelattice_out_%04d.csv", out_dir.c_str(), globalFrame);
            sprintf(meshfilename, "%s/DEM_particlelattice_%04d.vtk", out_dir.c_str(), globalFrame);

            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshfilename));

            DEMSim.DoDynamicsThenSync(frame_time); // Advance the simulation by frame_time

            globalFrame++; // Increment the global frame counter
        }
    }
    return 0;    
}