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

    const float R = 0.05;
    const float x_spacing = 2 * R;
    const float y_spacing = sqrt(3) * R;
    const int layers_z = 30;

    //prepare vectors for positions and a shared clump type
    std::vector<float3> positions;
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    //lattice generation
    for (int z = 0; z < layers_z; ++z) {
    int rows_xy = (z % 2 == 0) ? 61 : 60; // Alternating rows in the xy plane
    for (int y = 0; y < rows_xy; ++y) {
        for (int x = 0; x < rows_xy; ++x) {
            // Calculate position with an offset for every other row
            float xPos = x * x_spacing + ((y % 2) * R);
            float yPos = y * y_spacing;
            float zPos = z * 2 * R;
            positions.push_back(make_float3(xPos, yPos, zPos));
            clump_types.push_back(sph_type_1); // Assuming sph_type_1 is defined earlier as a sphere type
        }
    }}

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
    out_dir += "/lattice_validation";
    create_directory(out_dir);

    const int totalsteps = 5000;
    int outputFrequency = 100;

    for(int step = 0; step < totalsteps; ++step) {
        DEMSim.DoDynamicsThenSync(2e-5);
        if (step % outputFrequency == 0) {
        std::cout << "Simulation Step: " << step << std::endl;

        // Output the current state of the simulation
        char filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), step / outputFrequency);
        DEMSim.WriteSphereFile(std::string(filename));
    }
    
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowTimingStats();
    std::cout << "Simulation complete. Particles have settled" <<std::endl;
    return 0;
    }
}