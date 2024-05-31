//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This script tests the contact chain propagation within a DEM-medium, generated
// by an external force applied to the topmost element located on the plane symmetry
// for the system.
// written by btagliafierro@gmail.com May 23, 2024
// =============================================================================

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
    DEMSim.UseFrictionalHertzianModel();
    DEMSim.SetVerbosity("ERROR");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetOutputContent({"ABSV"});
    DEMSim.SetMeshOutputFormat("VTK");
    DEMSim.SetContactOutputContent(DEME_POINT | OWNER | FORCE | CNT_WILDCARD);

    std::cout << "============================================================" << std::endl;
    std::cout << "Initializing DEMdemo_ContactChain demo." << std::endl;
    
    float innerFriction = 0.7;
    float massMultiplier = 5.0;  // Magnitude of the external force
    
    std::cout << "Inner friction: " << innerFriction << "; Mass multiplier: " << massMultiplier << "." << std::endl;
    path out_dir = "";
    out_dir += "./ContactChain_impact_out";
    remove_all(out_dir);
    create_directories(out_dir);

    // E, nu, CoR, mu, Crr...
    auto mat_type_terrain =
        DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.33}, {"CoR", 0.3}, {"mu", innerFriction}, {"Crr", 0.0}});
    // this second material definition is used at time zero per the particle initialization.
    auto mat_type_terrain_zero =
        DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.33}, {"CoR", 0.0}, {"mu", innerFriction}, {"Crr", 0.0}});

    // Defining some of the quantities that are used later for this script.
    float terrain_rad = 0.01;
    float gravityMagnitude = 9.81;
    // See Equation (4) in Zhang et al. (2024).
    float DTc = PI * terrain_rad * sqrt(1e3 / (1e7 / (2 * (1 + 0.33)))) / (0.8766 + 0.163 * 0.33);
    // float step_size = 5 * 1.0e-5 * sqrt(terrain_rad * gravityMagnitude);
    float step_size = 0.1 * DTc;
    double world_sizeX = 122.0 * terrain_rad;  // 122.0 works fine for frictionless
    double world_sizeZ = 27 * terrain_rad;

    DEMSim.InstructBoxDomainDimension({-world_sizeX / 2., world_sizeX / 2.}, {-5 * terrain_rad, 5 * terrain_rad},
                                      {-1 * world_sizeZ, 10 * terrain_rad});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Creating the two clump templates we need, which are just spheres
    std::vector<std::shared_ptr<DEMClumpTemplate>> templates_terrain;

    templates_terrain.push_back(DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 4 / 3 * 1.0e3 * PI,
                                                      terrain_rad, mat_type_terrain_zero));

    templates_terrain.push_back(DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 4 / 3 * 1.0e3 * PI,
                                                      terrain_rad, mat_type_terrain));

    // Loading the position of the spheres from an external file.
    //! Note that this list does not include the particle located at (0.0,0.0).
    auto data_xyz = DEMSim.ReadClumpXyzFromCsv("./data/clumps/ContactChain_initial.csv");
    std::vector<float3> input_xyz;

    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;

    std::cout << data_xyz.size() << " Data points are loaded from the external list." << std::endl;

    for (unsigned int i = 0; i < data_xyz.size(); i++) {
        char t_name[20];
        sprintf(t_name, "%d", i);

        auto this_type_xyz = data_xyz[std::string(t_name)];
        input_xyz.insert(input_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());

        input_pile_template_type.push_back(templates_terrain[0]);
    }

    auto allParticles = DEMSim.AddClumps(input_pile_template_type, input_xyz);
    allParticles->SetFamily(1);

    // Separately, we include here the particle at (0.0,0.0), which is the one that will proxy the exterbal force.
    auto zeroParticle = DEMSim.AddClumps(templates_terrain[1], make_float3(0, 0, -terrain_rad));
    zeroParticle->SetFamily(3);
    auto driver = DEMSim.Track(zeroParticle);

    // Prescribe the acceleration before initializing
    float Aext = -gravityMagnitude * (massMultiplier);
    DEMSim.AddFamilyPrescribedAcc(2, "none", "none", std::to_string(Aext)); // Prescribed acceleration

    std::cout << "Total num of particles: " << (int)input_pile_template_type.size() + 1 << "." << std::endl;

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -gravityMagnitude));

    DEMSim.Initialize();

    float sim_time = 20.0;
    float time_settling = 5.0;
    unsigned int fps = 5;
    float frame_time = 1.0 / fps;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Selected output at " << fps << " FPS." << std::endl;
    unsigned int currframe = 0;
    double terrain_max_z;

    bool changeMaterial = true;
    bool impact_applied;

    for (float t = 0; t < sim_time; t += frame_time) {
        std::cout << "Output file: " << currframe << " at time " << t << " s." << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);

        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteContactFile(std::string(cnt_filename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);
        // Half time through the settling, the targeted properties for the material are applied.
        if (t > time_settling / 2 && changeMaterial) {
            DEMSim.DoDynamicsThenSync(0);
            std::cout << "Including restitution coefficient." << std::endl;
            changeMaterial = false;
            DEMSim.SetFamilyClumpMaterial(1, mat_type_terrain);
        }
        // Apply the impact load by changing the family
        if (t > time_settling && !impact_applied) {
            DEMSim.DoDynamicsThenSync(0);
            DEMSim.ChangeFamily(3, 2);
            impact_applied = true;
            std::cout << "Impact load applied at time " << t << " s." << std::endl;
        }
    }

    DEMSim.ShowTimingStats();
    std::cout << "==============================================================" << std::endl;
    std::cout << "DEMdemo_ContactChain exiting..." << std::endl;
    return 0;
}
