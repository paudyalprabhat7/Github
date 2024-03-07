//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

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

///////////////////////////////////////////////////////////////////////////////////////////////////////
void write_csv(std::string filename, std::vector<std::pair<std::string, std::vector<float>>> dataset){
    // Make a CSV file with one or more columns of integer values
    // Each column of data is represented by the pair <column name, column data>
    //   as std::pair<std::string, std::vector<int>>
    // The dataset is represented as a vector of these columns
    // Note that all columns should be the same size
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < dataset.size(); ++j)
    {
        myFile << dataset.at(j).first;
        if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << "\n";
    
    // Send data to the stream
    for(int i = 0; i < dataset.at(0).second.size(); ++i)
    {
        for(int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).second.at(i);
            if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }
    
    // Close the file
    myFile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// | OUTPUT_CONTENT:: VEL //|  OUTPUT_CONTENT:: ACC 
// | NORMAL


int main() {

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV  |  OUTPUT_CONTENT:: ABS_ACC );
    DEMSim.SetContactOutputContent(OWNER | FORCE | POINT | COMPONENT  );
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.EnsureKernelErrMsgLineNum();

    // E, nu, CoR, mu, Crr...
    auto mat_type_screw = DEMSim.LoadMaterial({{"E", 4.1e9}, {"nu", 0.394}, {"CoR", 0.3}, {"mu", 0.04}, {"Crr", 0.01}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 100e9}, {"nu", 0.45}, {"CoR", 0.5}, {"mu", 0.4}, {"Crr", 0.175}});
    auto mat_type_container = DEMSim.LoadMaterial({{"E", 100e9}, {"nu", 0.45}, {"CoR", 0.5}, {"mu", 0.4}, {"Crr", 0.175}});
    // If you don't have this line, then CoR,Crr,mu between screw material and granular material will be average of the two.
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_screw, mat_type_terrain, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_screw, mat_type_terrain, 0.01);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_screw, mat_type_terrain, 0.04);

    DEMSim.SetMaterialPropertyPair("CoR", mat_type_container, mat_type_terrain, 0.36);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_container, mat_type_terrain, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_container, mat_type_terrain, 1.00);
    
    // 1- Boundary:
    float step_size = 1e-6;
    //float Sim_CED = 0; //Pa

    double xBoundary = 0.32;
    double yBoundary = 0.32;
    double zBoundary = 0.5;
    double safety_margin= 0.05;

    DEMSim.InstructBoxDomainDimension(xBoundary, yBoundary, zBoundary + safety_margin);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_container);
    auto bot_wall = DEMSim.AddBCPlane(make_float3(0,0,-zBoundary/2.0),make_float3(0,0,1),mat_type_container);
    DEMSim.SetErrorOutVelocity(100.);

    
    //2-Bottom Wall:
    bot_wall->SetMass(1000);//1000 is the assumed mass of boundary. 
    bot_wall->SetMOI(make_float3(1./12.*1000.0*yBoundary*yBoundary, 1./12.*1000.0*xBoundary*xBoundary, 1./12.*1000.0*(xBoundary*xBoundary+yBoundary*yBoundary)));
    auto bot_wall_tracker = DEMSim.Track(bot_wall);

        // 2- Screw:

    // At the beginning of the simulation and until the end of the copression, the screw has a family type =1. this family is set to be fixed and does not contact the particles.
    // After, the screw family type is set to 3. 
    auto projectile = DEMSim.AddWavefrontMeshObject("./screwCM.obj", mat_type_screw);
    std::cout << "Total num of triangles: " << projectile->GetNumTriangles() << std::endl;
    projectile->InformCentroidPrincipal(make_float3(0.095965,0.237612,0.095964),make_float4 ( 0,0,0,1));
    projectile->SetInitPos(make_float3(0, 0.061, 0.1 ));
    float screw_mass = 0.59; // screw mass 0.59kg . 
    float I_XX = 0.0062;
    float I_YY = 0.003;
    float I_ZZ = 0.0062;
    projectile->SetMass(screw_mass);
    projectile->SetMOI(make_float3(I_XX, I_YY, I_ZZ));
    projectile->SetFamily(1);
    DEMSim.SetFamilyFixed(1);
    DEMSim.DisableContactBetweenFamilies(0, 1);
    // Track the projectile
    auto proj_tracker = DEMSim.Track(projectile);
    
    float rev_per_sec = 2;
    float ang_vel_Z = rev_per_sec * 3.14;
    float w_r = -3.0; //rad/s

    // Screw Coredinate system with respect to global corrdinate system:
    // y-aixs = y-axis global.
    // z-aixs = z-axis global.
    // x-aixs = x-axis global.
    float Pull_Force = 0; //N
    // In fact, because the cone's motion is completely pre-determined, we can just prescribe family 3
    DEMSim.SetFamilyPrescribedLinVel(2,"0", "0",  "none", false);// drop screw under gravity
    DEMSim.AddFamilyPrescribedAcc(2,  "none","none", to_string_with_precision(-100.0/screw_mass));
    DEMSim.SetFamilyPrescribedLinVel(3,  "0","none", "none",  false);
    DEMSim.AddFamilyPrescribedAcc(3,  "none","none", to_string_with_precision(-100.0/screw_mass));
    DEMSim.SetFamilyPrescribedAngVel(3,  "0", to_string_with_precision(w_r), "0",  false);


    //This code generates a normal distribution of clumps with a specified mean and standard deviation. 
    //Then the code samples a 1000 values of the distribution to generate clumps. 
    //The code in a way it ignores any sampled value that is less or more than the mean -/+ standered deviation. 
    // 3-Particles:
    float mean_radius = 0.0015;
    float std_radius = 0.0005;
    std::default_random_engine generator;
    std::normal_distribution<float> distribution(mean_radius, std_radius);//Normal Dist. of particles
    int num_particles = 4000; // number if different clumbs types (number of samples taken from distribution)
    //auto template_terrain = DEMSim.LoadSphereType(0.0, 0.0, mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    for (int i = 0; i < num_particles; i++) {
    float radius = distribution(generator);
    while(radius < mean_radius - std_radius || radius > mean_radius + std_radius) {
        radius = distribution(generator); // sample a different radius value until the condition is false
    }
    auto clump_template = DEMSim.LoadSphereType(radius * radius * radius * 2.5e3 * 4.0 / 3.0 * 3.14, radius, mat_type_terrain);
    clump_types.push_back(clump_template);
}
// Generate initial clumps for piling
    float spacing = 2.005*(mean_radius+std_radius);
    float sample_halfwidth_x = xBoundary / 2.000  - 2.005*(mean_radius+std_radius);
    float sample_halfwidth_y = yBoundary / 2.000  - 2.005*(mean_radius+std_radius);
    float fill_height = 0.5;
    float fill_bottom = -zBoundary / 2.0+ spacing;
    PDSampler sampler(spacing);
// Use a PDSampler-based clump generation process
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
    std::vector<float3> input_pile_xyz;
    float layer_z = 0;
    while (layer_z < fill_height) {
    float3 sample_center = make_float3(0, 0, fill_bottom + layer_z );
    auto layer_xyz = sampler.SampleBox(sample_center, make_float3(sample_halfwidth_x, sample_halfwidth_y, 0));
    unsigned int num_clumps = layer_xyz.size();
    // Select from available clump types
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_pile_template_type.push_back(clump_types.at(i % num_particles));
        }
    input_pile_xyz.insert(input_pile_xyz.end(), layer_xyz.begin(), layer_xyz.end());
    layer_z += spacing;
    }
// Calling AddClumps a to add clumps to the system
    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);
    the_pile->SetFamily(0);
    
    
    std::cout << "Terrain loaded: " <<  std::endl;
    size_t n_particles = input_pile_xyz.size();
    std::cout << "Added Number of clumps:" << n_particles << std::endl;
    // Create a inspector to find out stuff
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
    float max_z;
    float min_z;
    auto KE_finder = DEMSim.CreateInspector("clump_kinetic_energy");
    float KE;

    // A custom force model can be read in through a file and used by the simulation. Magic, right?
    //auto my_force_model = DEMSim.ReadContactForceModel("SJKR.cu");
    // This custom force model still uses contact history arrays, so let's define it
    //my_force_model->SetPerContactWildcards({"delta_tan_x", "delta_tan_y", "delta_tan_z", "delta_time"});
    //std::cout << "Force Module loaded:" << n_particles << std::endl;
    // 4- Compressor:
    // Now add a plane to compress the sample
    //auto compressor = DEMSim.AddExternalObject();
    //compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, 1), mat_type_container);
    //compressor->SetFamily(10);
    //DEMSim.DisableContactBetweenFamilies(0, 10);
    //DEMSim.DisableContactBetweenFamilies(1, 10);
    //DEMSim.DisableContactBetweenFamilies(2, 10);
    //DEMSim.DisableContactBetweenFamilies(3, 10);
    //DEMSim.SetFamilyFixed(10);
    //auto compressor_tracker = DEMSim.Track(compressor);


    // 5- Simulation initilization:
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEMSim.SetCDUpdateFreq(15);
    DEMSim.SetMaxVelocity(50.);
    DEMSim.SetExpandSafetyMultiplier(1.1);
    DEMSim.SetInitBinSize(4*mean_radius+std_radius);
    DEMSim.Initialize();


    // 6- Drop Particles:

    std::cout << "//////////////// Drop Particles //////////////////: " << std::endl; 

    path out_dir = current_path();
    out_dir += "/NormalDistBedSettlingWithScrew";
    create_directory(out_dir);

    float settling_end = 0.45;
    float drop_end = 0.85;
    //float compression_end = 0.75;
    float rolling_end = 3.5;

    unsigned int fps = 30;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float sim_time = 0.0;
    unsigned int total_frames = (unsigned int)((rolling_end )* fps);
    

    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    std::vector<float> ScrewXForceVector(total_frames);
    std::vector<float> ScrewYForceVector(total_frames);
    std::vector<float> ScrewZForceVector(total_frames);

    std::vector<float> ScrewXtorqueVector(total_frames);
    std::vector<float> ScrewYtorqueVector(total_frames);
    std::vector<float> ScrewZtorqueVector(total_frames);

    std::vector<float> ScrewXVelVector(total_frames);
    std::vector<float> ScrewYVelVector(total_frames);
    std::vector<float> ScrewZVelVector(total_frames);

    std::vector<float> ScrewXVector(total_frames);
    std::vector<float> ScrewYVector(total_frames);
    std::vector<float> ScrewZVector(total_frames);
    std::vector<float> timeVector(total_frames);

    std::vector<float> KEVector(total_frames);

    std::vector<float> BC_XForceVector(total_frames);
    std::vector<float> BC_YForceVector(total_frames);
    std::vector<float> BC_ZForceVector(total_frames);


    float3 SCREW_position;
    


    
    
    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    for (float t = 0; t < settling_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
        std::cout << "Frame: " << currframe<< " of " << total_frames << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        //sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
         
        
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        //DEMSim.WriteContactFile(std::string(cnt_filename));
        float3 pos_screw = proj_tracker->Pos();
        float3 force = proj_tracker->ContactAcc();
        float3 torque_screw = proj_tracker->ContactAngAccLocal();
        float3 VelocityScrew = proj_tracker->Vel();
        KE = KE_finder->GetValue();


        std::cout << "Max velocity of clump " << max_v_finder->GetValue() <<  " m/s" << std::endl;

        ScrewXForceVector[currframe] = force.x;
        ScrewYForceVector[currframe] = force.y;
        ScrewZForceVector[currframe] = force.z;

        ScrewXtorqueVector[currframe] = torque_screw.x*I_XX;
        ScrewYtorqueVector[currframe] = torque_screw.y*I_YY;
        ScrewZtorqueVector[currframe] = torque_screw.z*I_ZZ;

        ScrewXVelVector[currframe] = VelocityScrew.x;
        ScrewYVelVector[currframe] = VelocityScrew.y;
        ScrewZVelVector[currframe] = VelocityScrew.z;

        timeVector[currframe] = sim_time;
        ScrewXVector[currframe] = pos_screw.x;
        ScrewYVector[currframe] = pos_screw.y;
        ScrewZVector[currframe] = pos_screw.z;

        KEVector[currframe] = KE;

        float3 BC_pos = bot_wall_tracker->Pos();
        float3 BC_force = (bot_wall_tracker->ContactAcc())/1000.0;
        //std::cout << "Bottom wall pos: " << BC_pos.x << ", " << BC_pos.y << ", " << BC_pos.z << std::endl;
        //std::cout << "Bottom wall force: " << BC_force.x << ", " << BC_force.y << ", " << BC_force.z << std::endl;
        BC_XForceVector[currframe] = BC_force.x;
        BC_YForceVector[currframe] = BC_force.y;
        BC_ZForceVector[currframe] = BC_force.z;

        std::cout << "Time: " << sim_time << std::endl;
        std::cout << "Force on SCREW: " << force.x << ", " << force.y << ", " << force.z << std::endl;
       
        currframe++;

        }
        sim_time +=step_size;
        DEMSim.DoDynamics(step_size);   
    }
 

    max_z = max_z_finder->GetValue();
    min_z = min_z_finder->GetValue();
    std::cout << "Max Z is: " << max_z << std::endl;
    std::cout << "Min Z is: " << min_z << std::endl;
    

    // 7-Compresse:

    //std::cout << "//////////////// Compression //////////////////: " << std::endl;

    //DEMSim.EnableContactBetweenFamilies(0, 10);
    //double compressor_vel = 0.05;
    //float terrain_max_z = max_z_finder->GetValue();
    //compressor_tracker->SetPos(make_float3(0, 0, terrain_max_z));
    //DEMSim.DoDynamicsThenSync(step_size);

   //for (float t = sim_time; t < compression_end; t += step_size, curr_step++) {
    //if (curr_step % out_steps == 0) {
    //    std::cout << "Frame: " << currframe<< " of " << total_frames << std::endl;
    //    char filename[200], meshfilename[200], cnt_filename[200];
    //    sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
    //    sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
    //    //sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
    //    DEMSim.WriteSphereFile(std::string(filename));
    //    DEMSim.WriteMeshFile(std::string(meshfilename));
     //   float3 pos_screw = proj_tracker->Pos();
     //   float3 force = proj_tracker->ContactAcc();
     //   float3 torque_screw = proj_tracker->ContactAngAccLocal();
     //   float3 VelocityScrew = proj_tracker->Vel();
     //   KE = KE_finder->GetValue();


     //   std::cout << "Max velocity of clump " << max_v_finder->GetValue() <<  " m/s" << std::endl;

     //   ScrewXForceVector[currframe] = force.x;
     //   ScrewYForceVector[currframe] = force.y;
      //  ScrewZForceVector[currframe] = force.z;

     //   ScrewXtorqueVector[currframe] = torque_screw.x*I_XX;
     //   ScrewYtorqueVector[currframe] = torque_screw.y*I_YY;
     //   ScrewZtorqueVector[currframe] = torque_screw.z*I_ZZ;

     //   ScrewXVelVector[currframe] = VelocityScrew.x;
     //   ScrewYVelVector[currframe] = VelocityScrew.y;
     //   ScrewZVelVector[currframe] = VelocityScrew.z;

     //   timeVector[currframe] = sim_time;
     //   ScrewXVector[currframe] = pos_screw.x;
    //   ScrewYVector[currframe] = pos_screw.y;
    //    ScrewZVector[currframe] = pos_screw.z;

    //    KEVector[currframe] = KE;

    //    float3 BC_pos = bot_wall_tracker->Pos();
    //    float3 BC_force = (bot_wall_tracker->ContactAcc())/1000.0;
        //std::cout << "Bottom wall pos: " << BC_pos.x << ", " << BC_pos.y << ", " << BC_pos.z << std::endl;
        //std::cout << "Bottom wall force: " << BC_force.x << ", " << BC_force.y << ", " << BC_force.z << std::endl;
    //    BC_XForceVector[currframe] = BC_force.x;
    //    BC_YForceVector[currframe] = BC_force.y;
    //    BC_ZForceVector[currframe] = BC_force.z;

    //    std::cout << "Time: " << sim_time << std::endl;
        //std::cout << "Force on SCREW: " << force.x << ", " << force.y << ", " << force.z << std::endl;
        //std::cout << "Drawbar pull coeff: " << force.y / (screw_mass*9.81) << std::endl;

    //    currframe++;
    //}

    //    terrain_max_z -= compressor_vel * step_size;
    //    compressor_tracker->SetPos(make_float3(0, 0, terrain_max_z));
    //    sim_time +=step_size;
    //    DEMSim.DoDynamicsThenSync(step_size); 

    //}

    //Remove compressor
    
    //DEMSim.DisableContactBetweenFamilies(0, 10);
    //DEMSim.DisableContactBetweenFamilies(1, 10);
    //DEMSim.DisableContactBetweenFamilies(2, 10);
    //DEMSim.DisableContactBetweenFamilies(3, 10);
    //DEMSim.DoDynamicsThenSync(step_size);
    std::cout << " Simulation Time:  " << sim_time << std::endl;
    float terrain_max_z = max_z_finder->GetValue();
    proj_tracker->SetPos(make_float3(0, 0.061, terrain_max_z+0.1)); //0.2 screw diamerter
    DEMSim.DoDynamicsThenSync(step_size);
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();
    float3 initial_position = proj_tracker->Pos();
    max_z = max_z_finder->GetValue();
    min_z = min_z_finder->GetValue();
    std::cout << "Max Z is: " << max_z << std::endl;
    std::cout << "Min Z is: " << min_z << std::endl;
    std::cout << " Screw Postion Set at x = " << initial_position.x << " y = " << initial_position.y <<" z = " << initial_position.z << std::endl;

    // Drop Screw:   

    DEMSim.ChangeFamily(1, 2);
    DEMSim.DoDynamicsThenSync(step_size);


    std::cout << "//////////////// drop Screw //////////////////: " << std::endl; 

    for (float t = sim_time; t < drop_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
        std::cout << "Frame: " << currframe<< " of " << total_frames << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        //sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
        
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        //DEMSim.WriteContactFile(std::string(cnt_filename));
        
        float3 pos_screw = proj_tracker->Pos();
        float3 force = proj_tracker->ContactAcc();
        float3 torque_screw = proj_tracker->ContactAngAccLocal();
        float3 VelocityScrew = proj_tracker->Vel();
        force *= screw_mass;
        KE = KE_finder->GetValue();
        


        std::cout << "Max velocity of clump " << max_v_finder->GetValue() <<  " m/s" << std::endl;
        std::cout << "screw postion: " << pos_screw.x << ", " << pos_screw.y << ", " << pos_screw.z << std::endl;


        ScrewXForceVector[currframe] = force.x;
        ScrewYForceVector[currframe] = force.y;
        ScrewZForceVector[currframe] = force.z;

        ScrewXVelVector[currframe] = VelocityScrew.x;
        ScrewYVelVector[currframe] = VelocityScrew.y;
        ScrewZVelVector[currframe] = VelocityScrew.z;

        ScrewXtorqueVector[currframe] = torque_screw.x*I_XX;
        ScrewYtorqueVector[currframe] = torque_screw.y*I_YY;
        ScrewZtorqueVector[currframe] = torque_screw.z*I_ZZ;
        
        timeVector[currframe] = sim_time;
        ScrewXVector[currframe] = pos_screw.x;
        ScrewYVector[currframe] = pos_screw.y;
        ScrewZVector[currframe] = pos_screw.z;
        KEVector[currframe] = KE;

        float3 BC_pos = bot_wall_tracker->Pos();
        float3 BC_force = (bot_wall_tracker->ContactAcc())/1000.0;
        BC_XForceVector[currframe] = BC_force.x;
        BC_YForceVector[currframe] = BC_force.y;
        BC_ZForceVector[currframe] = BC_force.z;
        std::cout << "Time: " << sim_time << std::endl;
        currframe++;
        
        }
        sim_time +=step_size;
        
        DEMSim.DoDynamics(step_size);
        
        
    }

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    // Spen Screw:   

    DEMSim.ChangeFamily(2, 3);
    DEMSim.DoDynamicsThenSync(step_size);


    std::cout << "//////////////// Spen Screw //////////////////: " << std::endl; 

    for (float t = sim_time; t < rolling_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
        std::cout << "Frame: " << currframe<< " of " << total_frames << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        //sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
        
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        //DEMSim.WriteContactFile(std::string(cnt_filename));
        
        float3 pos_screw = proj_tracker->Pos();
        float3 force = proj_tracker->ContactAcc();
        float3 torque_screw = proj_tracker->ContactAngAccLocal();
        float3 VelocityScrew = proj_tracker->Vel();
        force *= screw_mass;
        KE = KE_finder->GetValue();
        
        
        std::cout << "Max velocity of clump " << max_v_finder->GetValue() <<  " m/s" << std::endl;
        std::cout << "screw postion: " << pos_screw.x << ", " << pos_screw.y << ", " << pos_screw.z << std::endl;


        ScrewXForceVector[currframe] = force.x;
        ScrewYForceVector[currframe] = force.y;
        ScrewZForceVector[currframe] = force.z;

        ScrewXVelVector[currframe] = VelocityScrew.x;
        ScrewYVelVector[currframe] = VelocityScrew.y;
        ScrewZVelVector[currframe] = VelocityScrew.z;

        ScrewXtorqueVector[currframe] = torque_screw.x*I_XX;
        ScrewYtorqueVector[currframe] = torque_screw.y*I_YY;
        ScrewZtorqueVector[currframe] = torque_screw.z*I_ZZ;
        
        timeVector[currframe] = sim_time;
        ScrewXVector[currframe] = pos_screw.x;
        ScrewYVector[currframe] = pos_screw.y;
        ScrewZVector[currframe] = pos_screw.z;
        float3 BC_pos = bot_wall_tracker->Pos();
        float3 BC_force = (bot_wall_tracker->ContactAcc())/1000.0;
        BC_XForceVector[currframe] = BC_force.x;
        BC_YForceVector[currframe] = BC_force.y;
        BC_ZForceVector[currframe] = BC_force.z;
        KEVector[currframe] = KE;
        std::cout << "Time: " << sim_time << std::endl;
        std::cout << "Force on SCREW: " << force.x << ", " << force.y << ", " << force.z << std::endl;
        

        currframe++;
        
        }
        sim_time += step_size;
        
        DEMSim.DoDynamics(step_size);
        
    }


    std::vector<std::pair<std::string, std::vector<float>>> vals = {{"Time", timeVector},{"PositionX", ScrewXVector},{"PositionY", ScrewYVector},{"PositionZ", ScrewZVector},{"Fx", ScrewXForceVector}, 
    {"Fy", ScrewYForceVector}, {"Fz", ScrewZForceVector},{"Vx", ScrewXVelVector}, {"Vy", ScrewYVelVector}, {"Vz", ScrewZVelVector},
    {"Tx", ScrewXtorqueVector}, {"Ty", ScrewYtorqueVector}, {"Tz", ScrewZtorqueVector},
    {"KE", KEVector}, {"BC_fx", BC_XForceVector}, {"BC_fy", BC_YForceVector}, {"BC_fz", BC_ZForceVector}};

   
    write_csv("Screw_Simulation_outputs_MixedP.csv",vals);

    
    std::cout << "simulation time: " << sim_time << std::endl;
    std::cout << "DEMdemo_ScrewDrop exiting..." << std::endl;
    
    
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();


    DEMSim.ShowTimingStats();
  
    return 0;

}
