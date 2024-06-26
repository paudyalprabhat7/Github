//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Enter code description here
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>

using namespace deme;
using namespace std::filesystem;
const std::string force_csv_header = "point_x,point_y,point_z,force_x,force_y,force_z";

//writing contactforces
void writeFloat3VectorsToCSV(const std::string& header,
                             const std::vector<std::vector<float3>>& vectors,
                             const std::string& filename,
                             size_t num_items);

int main() {
    DEMSolver DEMSim; //declare and initialize the object DEMSim of the DEMSolver class
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV); 
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.EnsureKernelErrMsgLineNum(); 

    //load the material properties
    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 2.1e11}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e11}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.65}, {"Crr", 0.01}});
    auto mat_type_analyticalb = DEMSim.LoadMaterial({{"E", 5e7}, {"nu", 0.5}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
    //auto mat_type_flexibleb = DEMSim.LoadMaterial({{"E", 3e6}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
    DEMSim.SetMaterialPropertyPair("mu", mat_type_terrain, mat_type_analyticalb, 0.5);

    //step size
    float terrain_rad = 0.01;
    float step_size = 1e-5;
    float world_size = 25 * terrain_rad;

    //Analytical boundary definition
    auto walls = DEMSim.AddExternalObject();
    auto bottom_wall = DEMSim.AddExternalObject();
    // Define each plane of the box domain manually
    bottom_wall->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, 1), mat_type_terrain); // Bottom plane
    //walls->AddPlane(make_float3(0, 0, world_size), make_float3(0, 0, -1), mat_type_analyticalb); // Top plane
    walls->AddPlane(make_float3(world_size / 2, 0, 0), make_float3(-1, 0, 0), mat_type_terrain); // Right plane
    walls->AddPlane(make_float3(-world_size / 2, 0, 0), make_float3(1, 0, 0), mat_type_terrain); // Left plane
    walls->AddPlane(make_float3(0, world_size / 2, 0), make_float3(0, -1, 0), mat_type_terrain); // Front plane
    walls->AddPlane(make_float3(0, -world_size / 2, 0), make_float3(0, 1, 0), mat_type_terrain); // Back plane

    auto bottom_tracker = DEMSim.Track(bottom_wall);

    float fill_height = 0.995 * world_size;
    

    //addition of impact cube
    float cube_thickness = 0.05 * world_size;
    float cube_size = 0.5 * world_size;
    float cube_pos = world_size + terrain_rad * 10;
    
    auto projectile = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_cube);
    projectile->Scale(make_float3(cube_size, cube_size, cube_thickness));

    projectile->SetInitPos(make_float3(0.0, 0.0, cube_pos));
    float cube_density = 7.6e3;
    float cube_mass = cube_density * (cube_thickness * cube_thickness * cube_size);
    projectile->SetMass(cube_mass);
    projectile->SetMOI(make_float3(cube_mass * 1 / 6, cube_mass * 1 / 6, cube_mass * 1 / 6));
    projectile->SetFamily(2);
    DEMSim.SetFamilyFixed(2);

    //terrain definition
    
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.69e3 * 4/3 * 3.141, terrain_rad, mat_type_terrain);

    //terrain sampling
    HCPSampler sampler(terrain_rad * 2.2);
    float3 fill_center = make_float3(0, 0, fill_height/2 + 2 * terrain_rad);
    float3 fill_halfsize = make_float3(world_size/2, world_size/2, fill_height/2);
    auto input_xyz = sampler.SampleBox(fill_center, fill_halfsize);
    auto particles = DEMSim.AddClumps(template_terrain, input_xyz);

    std::cout << "Total num of particles: " << particles->GetNumClumps() << std::endl;
    std::cout << "Total num of spheres: " << particles->GetNumSpheres() << std::endl;

    //initialization of simulation
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0 , -9.81));
    DEMSim.SetMaxVelocity(15.);

    DEMSim.Initialize();

    //creating the output directory
    path out_dir = current_path();
    out_dir += "/Output_rubber_10x";
    create_directory(out_dir);

    //visualization frame time
    float sim_time = 4.0;
    float settle_time = 2.0;
    unsigned int fps = 24;
    float frame_time = 1.0 / fps;
    std::cout << "Output at " << fps <<" FPS" << std::endl;

    std::vector<float3> forces, points;
    size_t num_force_pairs = 0;
    
    unsigned int curr_frame = 0;

    //loop for settling
    for (float t = 0; t<settle_time; t+=frame_time) {
        char filename[200], force_filename[200], meshfilename[200];
        sprintf(filename, "%s/simulation_output_%04d.csv", out_dir.c_str(), curr_frame);
        sprintf(meshfilename, "%s/simulation_output_%04d.vtk", out_dir.c_str(), curr_frame);
        sprintf(force_filename, "%s/DEMdemo_forces_%04d.csv", out_dir.c_str(), curr_frame);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        writeFloat3VectorsToCSV(force_csv_header, {points, forces}, force_filename, num_force_pairs);
        curr_frame++;
        num_force_pairs = bottom_tracker -> GetContactForces(points, forces);
        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();
    }

    //dropping the cube
    DEMSim.ChangeFamily(2,1);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (float t = 0; t < sim_time; t+=frame_time) {
        char filename[200], force_filename[200], meshfilename[200];
        sprintf(filename, "%s/simulation_output_%04d.csv", out_dir.c_str(), curr_frame);
        sprintf(meshfilename, "%s/simulation_output_%04d.vtk", out_dir.c_str(), curr_frame);
        sprintf(force_filename, "%s/DEMdemo_forces_%04d.csv", out_dir.c_str(), curr_frame);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        writeFloat3VectorsToCSV(force_csv_header, {points, forces}, force_filename, num_force_pairs);
        curr_frame++;
        num_force_pairs = bottom_tracker -> GetContactForces(points, forces);
        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();

    }

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    //post simulation housekeeping
    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();
    std::cout << "Simulation exiting" << std::endl;
    return 0;
}

void writeFloat3VectorsToCSV(const std::string& header,
                             const std::vector<std::vector<float3>>& vectors,
                             const std::string& filename,
                             size_t num_items) {
    std::ofstream file(filename);

    // Check if the file was successfully opened
    if (!file.is_open()) {
        std::cout << "Failed to open the CSV file used to store mesh forces!" << std::endl;
        return;
    }

    file << header << "\n";

    // Write vectors as columns
    for (size_t i = 0; i < num_items; ++i) {
        for (size_t j = 0; j < vectors.size(); ++j) {
            if (i < vectors[j].size()) {
                file << vectors[j][i].x << "," << vectors[j][i].y << "," << vectors[j][i].z;
            }
            if (j != vectors.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    // Close the file
    file.close();
}