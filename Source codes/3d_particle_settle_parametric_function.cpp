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
#include <string>

using namespace deme;
using namespace std::filesystem;

// CSV header for contact forces
const std::string force_csv_header = "point_x,point_y,point_z,force_x,force_y,force_z";

// Function to write contact forces to a CSV file
void writeFloat3VectorsToCSV(const std::string& header,
                             const std::vector<std::vector<float3>>& vectors,
                             const std::string& filename,
                             size_t num_items);

int main() {
    // Define parameter values for different simulations
    float bottom_boundary_E[] = {1e8, 2e8, 3e8, 4e8, 5e8};  // Elastic moduli for bottom boundary
    float side_planes_E[] = {1e7, 2e7, 3e7};  // Elastic moduli for side planes
    float drop_heights[] = {2.1, 2.2, 2.3};  // Drop heights for the cube

    // Create the master directory for all simulation results
    path master_dir = current_path() / "SimulationResults_trial15May2024";
    create_directories(master_dir);

    // Calculate the number of different parameter combinations
    size_t num_bottom_E = sizeof(bottom_boundary_E) / sizeof(bottom_boundary_E[0]);
    size_t num_side_E = sizeof(side_planes_E) / sizeof(side_planes_E[0]);
    size_t num_drop_heights = sizeof(drop_heights) / sizeof(drop_heights[0]);

    // Loop through all parameter combinations
    for (size_t i = 0; i < num_bottom_E; ++i) {
        for (size_t j = 0; j < num_side_E; ++j) {
            for (size_t k = 0; k < num_drop_heights; ++k) {
                // Encapsulate the simulation in a separate scope
                try {
                    // Get current parameter values
                    float E_bottom = bottom_boundary_E[i];
                    float E_side = side_planes_E[j];
                    float drop_height = drop_heights[k];

                    std::cout << "Starting simulation with E_bottom: " << E_bottom << ", E_side: " << E_side << ", drop_height: " << drop_height << std::endl;
                    
                    // Initialize the DEM solver
                    DEMSolver DEMSim; // Declare and initialize the DEMSolver object
                    DEMSim.SetVerbosity(INFO);
                    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV); 
                    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
                    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
                    DEMSim.EnsureKernelErrMsgLineNum(); 

                    // Load material properties with updated elastic moduli
                    auto mat_type_cube = DEMSim.LoadMaterial({{"E", 2.1e10}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
                    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.65}, {"Crr", 0.01}});
                    auto mat_type_analyticalb = DEMSim.LoadMaterial({{"E", E_bottom}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
                    auto mat_type_flexibleb = DEMSim.LoadMaterial({{"E", E_side}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
                    DEMSim.SetMaterialPropertyPair("mu", mat_type_terrain, mat_type_analyticalb, 0.5);

                    // Define simulation parameters
                    float step_size = 1e-4;  // Time step size
                    float world_size = 2;  // Size of the simulation world

                    // Define the analytical boundaries (box domain)
                    auto walls = DEMSim.AddExternalObject();
                    auto bottom_wall = DEMSim.AddExternalObject();
                    bottom_wall->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, 1), mat_type_analyticalb); // Bottom plane
                    walls->AddPlane(make_float3(world_size / 2, 0, 0), make_float3(-1, 0, 0), mat_type_flexibleb); // Right plane
                    walls->AddPlane(make_float3(-world_size / 2, 0, 0), make_float3(1, 0, 0), mat_type_flexibleb); // Left plane
                    walls->AddPlane(make_float3(0, world_size / 2, 0), make_float3(0, -1, 0), mat_type_flexibleb); // Front plane
                    walls->AddPlane(make_float3(0, -world_size / 2, 0), make_float3(0, 1, 0), mat_type_flexibleb); // Back plane

                    // Track contact forces on the bottom wall
                    auto bottom_tracker = DEMSim.Track(bottom_wall);

                    // Define the terrain fill height
                    float fill_height = 0.995 * world_size;

                    // Define the impact cube properties
                    float cube_thickness = 0.05 * world_size;
                    float cube_size = 0.5 * world_size;

                    // Add the impact cube to the simulation
                    auto projectile = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_cube);
                    projectile->Scale(make_float3(cube_size, cube_size, cube_thickness));
                    projectile->SetInitPos(make_float3(0.0, 0.0, drop_height));
                    float cube_density = 7.6e3;
                    float cube_mass = cube_density * (cube_size * cube_size * cube_thickness);
                    projectile->SetMass(cube_mass);
                    projectile->SetMOI(make_float3(cube_mass * 1 / 6, cube_mass * 1 / 6, cube_mass * 1 / 6));
                    projectile->SetFamily(2);
                    DEMSim.SetFamilyFixed(2);

                    // Define the terrain particles
                    float terrain_rad = 0.08;
                    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.69e3 * 4/3 * 3.141, terrain_rad, mat_type_terrain);

                    // Sample the terrain
                    HCPSampler sampler(terrain_rad * 2.2);
                    float3 fill_center = make_float3(0, 0, fill_height / 2 + 2 * terrain_rad);
                    float3 fill_halfsize = make_float3(world_size / 2, world_size / 2, fill_height / 2);
                    auto input_xyz = sampler.SampleBox(fill_center, fill_halfsize);
                    auto particles = DEMSim.AddClumps(template_terrain, input_xyz);

                    std::cout << "Total num of particles: " << particles->GetNumClumps() << std::endl;
                    std::cout << "Total num of spheres: " << particles->GetNumSpheres() << std::endl;

                    // Initialize the simulation
                    DEMSim.SetInitTimeStep(step_size);
                    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
                    DEMSim.SetMaxVelocity(15.);
                    DEMSim.Initialize();

                    // Create output directory based on parameters within the master directory
                    path out_dir = master_dir / ("BottomBoundary_E_" + std::to_string(static_cast<int>(E_bottom))) /
                                   ("SidePlanes_E_" + std::to_string(static_cast<int>(E_side))) /
                                   ("DropHeight_" + std::to_string(static_cast<int>(drop_height * 10)));
                    create_directories(out_dir);

                    // Simulation settings
                    float sim_time = 4.0;  // Simulation duration
                    float settle_time = 2.0;  // Settling time
                    unsigned int fps = 24;  // Frames per second for output
                    float frame_time = 1.0 / fps;
                    std::cout << "Output at " << fps << " FPS" << std::endl;

                    std::vector<float3> forces, points;
                    size_t num_force_pairs = 0;
                    unsigned int curr_frame = 0;

                    // Loop for settling
                    for (float t = 0; t < settle_time; t += frame_time) {
                        char filename[200], force_filename[200], meshfilename[200];
                        sprintf(filename, "%s/simulation_output_%04d.csv", out_dir.c_str(), curr_frame);
                        sprintf(meshfilename, "%s/simulation_output_%04d.vtk", out_dir.c_str(), curr_frame);
                        sprintf(force_filename, "%s/DEMdemo_forces_%04d.csv", out_dir.c_str(), curr_frame);
                        DEMSim.WriteSphereFile(std::string(filename));
                        DEMSim.WriteMeshFile(std::string(meshfilename));
                        writeFloat3VectorsToCSV(force_csv_header, {points, forces}, force_filename, num_force_pairs);
                        curr_frame++;
                        num_force_pairs = bottom_tracker->GetContactForces(points, forces);
                        DEMSim.DoDynamicsThenSync(frame_time);
                        DEMSim.ShowThreadCollaborationStats();
                    }

                    // Drop the cube
                    DEMSim.ChangeFamily(2, 1);

                    // Start timing the simulation
                    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
                    for (float t = 0; t < sim_time; t += frame_time) {
                        char filename[200], force_filename[200], meshfilename[200];
                        sprintf(filename, "%s/simulation_output_%04d.csv", out_dir.c_str(), curr_frame);
                        sprintf(meshfilename, "%s/simulation_output_%04d.vtk", out_dir.c_str(), curr_frame);
                        sprintf(force_filename, "%s/DEMdemo_forces_%04d.csv", out_dir.c_str(), curr_frame);
                        DEMSim.WriteSphereFile(std::string(filename));
                        DEMSim.WriteMeshFile(std::string(meshfilename));
                        writeFloat3VectorsToCSV(force_csv_header, {points, forces}, force_filename, num_force_pairs);
                        curr_frame++;
                        num_force_pairs = bottom_tracker->GetContactForces(points, forces);
                        DEMSim.DoDynamicsThenSync(frame_time);
                        DEMSim.ShowThreadCollaborationStats();
                    }

                    // End timing the simulation
                    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
                    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

                    // Post-simulation housekeeping
                    DEMSim.ShowTimingStats();
                    DEMSim.ShowAnomalies();
                    std::cout << "Simulation exiting" << std::endl;

                    // Clear vectors to free memory
                    forces.clear();
                    points.clear();
                    forces.shrink_to_fit();
                    points.shrink_to_fit();

                } catch (const std::bad_alloc& e) {
                    std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                    return 1;  // Exit if memory allocation fails
                } catch (const std::exception& e) {
                    std::cerr << "An error occurred: " << e.what() << std::endl;
                    return 1;  // Exit if any other exception occurs
                }
            }
        }
    }
    return 0;
}

// Function to write float3 vectors to a CSV file
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
