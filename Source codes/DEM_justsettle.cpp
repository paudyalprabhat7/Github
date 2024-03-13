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
    // Output less info at initialization
    DEMSim.SetVerbosity("ERROR");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetMeshOutputFormat("VTK");

    path out_dir = current_path();
    out_dir += "/DemoOutput_ParticleSettle";
    create_directory(out_dir);

    // Material properties for the terrain
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 7e7}, {"nu", 0.24}, {"CoR", 0.9}, {"mu", 0.3}, {"Crr", 0.0}});

    double terrain_rad = 0.006 / 2.;

    float step_size = 2e-6;
    double world_size = terrain_rad * 122.5;
    DEMSim.InstructBoxDomainDimension({-world_size / 2., world_size / 2.}, {-world_size / 2., world_size / 2.},
                                      {0, 10 * world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Generate terrain particles
    std::vector<std::shared_ptr<DEMClumpTemplate>> templates_terrain;
    for (int i = 0; i < 11; i++) {
        templates_terrain.push_back(DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.0e3 * 4 / 3 * PI,
                                                          terrain_rad, mat_type_terrain));
        terrain_rad += 0.0001 / 2.;
    }

    PDSampler sampler(2.01 * terrain_rad);
    float sample_halfwidth = world_size / 2 - 2 * terrain_rad;
    float fullheight = world_size * 6.;
    auto sample_center = make_float3(0, 0, fullheight / 2 + 1 * terrain_rad);
    auto input_xyz = sampler.SampleBox(sample_center, make_float3(sample_halfwidth, 0.f, fullheight / 2.));

    // Random selection of templates for each particle
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, templates_terrain.size() - 1);

    std::vector<std::shared_ptr<DEMClumpTemplate>> template_to_use(input_xyz.size());
    for (unsigned int i = 0; i < input_xyz.size(); i++) {
        template_to_use[i] = templates_terrain[dist(gen)];
    }
    DEMSim.AddClumps(template_to_use, input_xyz);

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

    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;
    std::cout << "Particle settling simulation completed." << std::endl;
    return 0;
}
