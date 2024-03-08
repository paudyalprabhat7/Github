//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// A meshed CUBE hitting a granular bed under gravity.
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

//lattice generation function check
std::vector<float3> GenerateCustomLattice(int n, int y, float particleDiameter) {
    std::vector<float3> positions;
    bool isNLayer = true; //n particles in the 1st layer
    float spacing = particleDiameter;

    for (int layer = 0; layer < y; ++layer) {
        int particlesInLayer = isNLayer ? n : n + 1;

        float startX = -particlesInLayer * spacing / 2.0f; 
        float startY = layer * spacing; 

        for (int i = 0; i < particlesInLayer; ++i) {
            positions.push_back({startX + i * spacing, startY, 0.0f}); 
        }
        isNLayer = !isNLayer; // Alternating between n and n+1 particles per layer
    }
    return positions;
}

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetErrorOutVelocity(20000.);
    DEMSim.EnsureKernelErrMsgLineNum();

    //material terrain
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 50e6}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.3}, {"Crr", 0.01}});
    
    //step size
    float step_size = 1e-6;

    //terrain radius and diameter
    float terrain_rad = 0.05;
    float particleDiameter = 2.0f * terrain_rad;

    //world size
    double world_size = 60 * terrain_rad;
    DEMSim.InstructBoxDomainDimension({0, world_size}, {0, world_size}, {0, world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.6e3 * 4 / 3 * 3.14,
                                                  terrain_rad, mat_type_terrain);

    
    //domaining the particle generation
    int n = 30; //base number of particles in x and y directions
    int y  = 2; //number of layers
    
    auto positions = GenerateCustomLattice(n, y, particleDiameter);

    //assigning particles to positions
    for(const auto&pos : positions) {
        DEMSim.AddClumps(template_terrain, pos);
    }

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
    // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
    // happen anyway and if it does, something already went wrong.
    DEMSim.SetMaxVelocity(15.);
    // In general you don't have to worry about SetExpandSafetyAdder, unless if an entity has the property that a point
    // on it can move much faster than its CoM. In this demo, you are dealing with a meshed ball and you in fact don't
    // have this problem. In the Centrifuge demo though, this can be a problem since the centrifuge's CoM is not moving,
    // but its pointwise velocity can be high, so it needs to be accounted for using this method.
    DEMSim.SetExpandSafetyAdder(5.);
    // You usually don't have to worry about initial bin size. In very rare cases, init bin size is so bad that auto bin
    // size adaption is effectless, and you should notice in that case kT runs extremely slow. Then in that case setting
    // init bin size may save the simulation.
    // DEMSim.SetInitBinSize(4 * terrain_rad);//

    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/DemoLayers";
    create_directory(out_dir);

    float settle_time = 2.0;
    unsigned int fps = 20;
    float frame_time = 1.0 / fps;

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;

    for (float t = 0; t < settle_time; t += frame_time) {
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();
    std::cout << "Simulation exiting" << std::endl;
    return 0;
}
