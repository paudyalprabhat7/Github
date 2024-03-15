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
    const int baseLayerCountX = 16; // specifying the number of particles in the bottom layer x-direction
    const int totalLayers = 2; // Total number of layers to stack

    std::vector<float3> positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Generate the custom lattice
    for (int layer = 0; layer < totalLayers; ++layer) {
        int layerCountX = baseLayerCountX - layer; // Decrease the number of particles for every other layer for a triangular lattice structure
        int layerCountY = static_cast<int>(std::round(0.1 * layerCountX)); // number of particles in the y-direction would essentially be 10% of that in the y-direction, rounded to the nearest integer
        for (int y = 0; y < layerCountY; ++y) {
            for (int x = 0; x < layerCountX; ++x) {
                
                //positions of particles in each layer in z-direction
                
                float xPos = x * particleDiameter;
                float yPos = y * sqrt(3) * R; // Keeping the close packing in y
                float zPos = layer * particleDiameter; // Stack the layers directly on top of each other

                positions.push_back(make_float3(xPos, yPos, zPos));
                clump_types.push_back(sph_type_1); // Assume sph_type_1 is defined as before
            }
        }
    }

    auto particles = DEMSim.AddClumps(clump_types, positions);
    
    // Add bottom plane mesh
    auto bot_plane = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/plane_20by20.obj").string(), mat_type_2);
    bot_plane->SetInitPos(make_float3(0, 0, -1.25));
    bot_plane->SetMass(10000.);

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