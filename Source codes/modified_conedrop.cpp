//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Validation study 1
// =============================================================================

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;

const double math_PI = 3.14159;

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent(OWNER | FORCE | POINT);
    DEMSim.SetErrorOutVelocity(20000.0);

    //Material and material pair definition
    // E, nu, CoR, mu, Crr...
    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.7}, {"Crr", 0.01}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.4}, {"Crr", 0.01}});
    // If you don't have this line, then values will take average between 2 materials, when they are in contact
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_cube, mat_type_terrain, 0.8);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_cube, mat_type_terrain, 0.7);

    float cube_speed = 0.03;
    float step_size = 5e-6;
    double world_size = 10;

    //adding world boundaries
    DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
    DEMSim.InstructBoxDomainBoundingBC("none", mat_type_terrain);
    
    //geometric parameters:
    float terrain_rad = 0.05; //radius of particles
    
    const float binWidth = 122 * terrain_rad; //bin width
    const float binLength = 122 * terrain_rad; //bin length 
    const float binHeight = 30 * terrain_rad; //bin height

    auto walls = DEMSim.AddExternalObject();

    double bottom = 0;
    
    /*
    // Bottom plane
    walls->AddPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);

    // Adding walls
    // Left wall
    walls->AddPlane(make_float3(-binLength / 2, 0, bottom + binHeight / 2), make_float3(1, 0, 0), mat_type_terrain);
    // Right wall
    walls->AddPlane(make_float3(binLength / 2, 0, bottom + binHeight / 2), make_float3(-1, 0, 0), mat_type_terrain);
    // Front wall
    walls->AddPlane(make_float3(0, -binWidth / 2, bottom + binHeight / 2), make_float3(0, 1, 0), mat_type_terrain);
    // Back wall
    walls->AddPlane(make_float3(0, binWidth / 2, bottom + binHeight / 2), make_float3(0, -1, 0), mat_type_terrain);
    */

    // adding the loading plate
    float dropobj_thickness = terrain_rad * 4.0;
    float dropobj_width = binWidth * 0.95;
    cube_speed = 0.03;

    auto projectile = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_cube);
    std::cout << "Total num of triangles: " << projectile->GetNumTriangles() << std::endl;

    projectile->Scale(make_float3(dropobj_width, dropobj_width, dropobj_thickness));

    float cube_mass = 7.8e3;
    projectile->SetMass(cube_mass);
    projectile->SetMOI(make_float3(cube_mass * 1 / 6, cube_mass * 1 / 6, cube_mass * 1 / 6));
    projectile->SetFamily(2);
    DEMSim.SetFamilyFixed(2);

    // Define the dimensions of the area to be filled
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.6e3 * 4 / 3 * 3.14,
                                                  terrain_rad, mat_type_terrain);

    float sample_halfwidth = binWidth / 2;
    float sample_halflength = binLength / 2;
    float sample_halfheight = binHeight / 2;

    // Center position for the terrain fill (assuming the bottom of the bin is at z=0)
    float3 sample_center = make_float3(0, 0, sample_halfheight);

    // Generate positions using cubic packing
    auto input_xyz = DEMBoxHCPSampler(sample_center, make_float3(sample_halfwidth, sample_halflength, sample_halfheight), 2.01 * terrain_rad);
    DEMSim.AddClumps(template_terrain, input_xyz);
    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;

    auto proj_tracker = DEMSim.Track(projectile);

    // Because the cube's motion is completely pre-determined, we can just prescribe family 1
    DEMSim.SetFamilyPrescribedLinVel(1, "0", "0", "-" + to_string_with_precision(cube_speed));
    // Initially cube is inactive and is assigned a family of 1
    DEMSim.SetFamilyFixed(2);
    DEMSim.DisableContactBetweenFamilies(0, 2);

    // Now add a plane to compress the sample
    auto compressor = DEMSim.AddExternalObject();
    compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
    compressor->SetFamily(10);
    DEMSim.SetFamilyFixed(10);
    auto compressor_tracker = DEMSim.Track(compressor);

    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    // auto total_volume_finder = DEMSim.CreateInspector("clump_volume", "return (X * X + Y * Y <= 0.25 * 0.25) && (Z <=
    // -0.3);");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // CD freq will be auto-adapted so it does not matter much here.
    DEMSim.SetCDUpdateFreq(20);
    // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
    // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
    // happen anyway and if it does, something already went wrong.
    DEMSim.SetMaxVelocity(10.);

    DEMSim.Initialize();

    std::filesystem::path out_dir = std::filesystem::current_path();
    out_dir += "/validation1_output";
    std::filesystem::create_directory(out_dir);

    // Settle
    DEMSim.DoDynamicsThenSync(0.8);

    // Compress until dense enough
    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    unsigned int fps = 20;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    double compressor_vel = 0.05;
    float terrain_max_z = max_z_finder->GetValue();
    double init_max_z = terrain_max_z;
    float bulk_density = -10000.;

    while (bulk_density < 1500.) {
        float matter_mass = total_mass_finder->GetValue();
        // Adjusted volume calculation for rectangular domain
        float total_volume = binWidth * binLength * (terrain_max_z - bottom);
        bulk_density = matter_mass / total_volume;

        if (curr_step % out_steps == 0) {
            char filename[200], meshname[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            // DEMSim.WriteMeshFile(std::string(meshname));
            std::cout << "Compression bulk density: " << bulk_density << std::endl;
            currframe++;
        }

        terrain_max_z -= compressor_vel * step_size;
        compressor_tracker->SetPos(make_float3(0, 0, terrain_max_z));
        DEMSim.DoDynamics(step_size);
        curr_step++;
    }
    // Then gradually remove the compressor
    while (terrain_max_z < init_max_z) {
        if (curr_step % out_steps == 0) {
            char filename[200], meshname[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            // DEMSim.WriteMeshFile(std::string(meshname));
            float matter_mass = total_mass_finder->GetValue();
            float total_volume = binWidth * binLength * (terrain_max_z - bottom);
            bulk_density = matter_mass / total_volume;
            std::cout << "Compression bulk density: " << bulk_density << std::endl;
            currframe++;
        }

        terrain_max_z += compressor_vel * step_size;
        compressor_tracker->SetPos(make_float3(0, 0, terrain_max_z));
        DEMSim.DoDynamics(step_size);
        curr_step++;
    }

    // Remove compressor
    DEMSim.DoDynamicsThenSync(0.);
    DEMSim.DisableContactBetweenFamilies(0, 10);
    DEMSim.DoDynamicsThenSync(0.2);
    terrain_max_z = max_z_finder->GetValue();

    float sim_end = 7.0;
    fps = 2500;
    float frame_time = 1.0 / fps;
    std::cout << "Output at " << fps << " FPS" << std::endl;

    // Put the cube in place
    double starting_height = terrain_max_z + terrain_rad*10;
    // Its initial position should be right above the cone tip...
    proj_tracker->SetPos(make_float3(0, 0, starting_height));

    // Enable cube
    DEMSim.ChangeFamily(2, 1);
    float matter_mass = total_mass_finder->GetValue();
    float total_volume = binWidth * binLength * (terrain_max_z - bottom);
    std::cout << "Bulk density: " << bulk_density << std::endl;

    bool contact_made = false;
    unsigned int frame_count = 0;
    double plate_z_when_first_contact;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (float t = 0; t < sim_end; t += frame_time) {
        //tracking the pressing plate if necessary
        float3 forces = proj_tracker->ContactAcc(); 
        float pressure = std::abs(forces.z) / (dropobj_width * dropobj_width); 

        // Condition to check for initial contact, if not already done
        if (!contact_made && pressure > 1e-4) { // Threshold might need adjustment
            contact_made = true;
            plate_z_when_first_contact = proj_tracker->GetPos()[2]; // Assuming proj_tracker provides position
        }

        // Calculate penetration (or compression depth) differently for pressing action
        float penetration = contact_made ? (plate_z_when_first_contact - proj_tracker->GetPos()[2]) : 0.0f;

        std::cout << "Time: " << t << std::endl;
        std::cout << "Plate Z coord: " << proj_tracker->GetPos()[2] << std::endl;
        std::cout << "Compression: " << penetration << std::endl;
        std::cout << "Force on plate: " << forces.x << ", " << forces.y << ", " << forces.z << std::endl;
        std::cout << "Pressure: " << pressure << std::endl;

        if (frame_count % 500 == 0) {
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));
            DEMSim.ShowThreadCollaborationStats();
        }

        DEMSim.DoDynamicsThenSync(frame_time);

        // Update for pressing plate's movement if it's being controlled dynamically
        frame_count++;
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    std::cout << "Simulation exiting..." << std::endl;
    return 0;
}