/**
 * @file 	hydrostatic_fsi.cpp
 * @brief 	structure deformation due to hydrostatic pressure under gravity.
 * @details This is the one of the basic test cases
 * for understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Yujie Zhu, Chi Zhang and Xiangyu Hu
 * @version 0.1
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI (3.14159265358979323846)
/*@brief Basic geometry parameters and numerical setup.
 */

Real resolution_ref = 0.008; /* Initial particle spacing*/

/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.2, -0.05, -0.2), Vec3d(0.2, 1.0, 0.2));

/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;                   /*Fluid density*/
Real rho0_a = 1.226;                    /*Air density*/
Real gravity_g = 9.81;                  /*Gravity force of fluid*/
Real U_f = 2.0 * sqrt(gravity_g * 0.5); /**< Characteristic velocity. */
Real U_g = 2.0 * sqrt(gravity_g * 0.5); /**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);       /**< Reference sound speed. */
Real f = 1.7;
Real a = 0.04;
Real c_p_water = 4.179e3;
Real c_p_air = 1.012e3;
Real k_water = 0.620;
Real k_air = 0.0254;
Real diffusion_coff_water = k_water / (c_p_water * rho0_f);
Real diffusion_coff_air = k_air / (c_p_air * rho0_a);
Real mu_water = 653.9e-6;
Real mu_air = 20.88e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.0, 0);

std::string fuel_tank_outer = "./input/tank_outer.STL";
std::string fuel_tank_inner = "./input/tank_inner.STL";
std::string water_05 = "./input/water_05.STL";
std::string air_05 = "./input/gas_05.STL";
std::string probe_s1_shape = "./input/ProbeS1.STL";
std::string probe_s2_shape = "./input/ProbeS2.STL";
std::string probe_s3_shape = "./input/ProbeS3.STL";

//std::string fuel_tank_inner = "./input/tank_inner_2_hole.STL";
//std::string water_05 = "./input/water_2_hole.STL";
/*
Fuel Tank.
*/
class Tank : public ComplexShape
{
  public:
    explicit Tank(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */

        add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale, "OuterWall");
        subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale, "InnerWall");
    }
};
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
    }
};

class VariableGravity : public Gravity
{
    Real time_ = 0;

  public:
    VariableGravity(Vecd gravity_vector) : Gravity(gravity_vector) {};
    virtual Vecd InducedAcceleration(const Vecd &position = Vecd::Zero()) override
    {
        time_ = GlobalStaticVariables::physical_time_;
        if (time_ <= 2.0)
        {
            global_acceleration_ = global_acceleration_;
        }
        else
        {
            global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
        }

        return global_acceleration_;
    }
};

class ProbeS1 : public ComplexShape
{
  public:
    explicit ProbeS1(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_probe(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
    }
};

class ProbeS2 : public ComplexShape
{
  public:
    explicit ProbeS2(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_probe_2(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
    }
};

class ProbeS3 : public ComplexShape
{
  public:
    explicit ProbeS3(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_probe_3(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
    }
};

StdVec<Vecd> liquid_temperature_observation_location = {Vecd(0.0, 0.35, 0.0), Vecd(0.0, 0.5, 0.0), Vecd(0.0, 0.6, 0.0)};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(false);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(true);
    /* Tag for computation from restart files. 0: start with initial condition. */
    sph_system.setRestartStep(0);
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment in_output(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody tank(sph_system, makeShared<Tank>("Tank"));
    //tank.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    tank.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? tank.generateParticles<BaseParticles, Reload>(tank.getName())
        : tank.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of gate observer.
    //----------------------------------------------------------------------
    ObserverBody liquid_temperature_observer(sph_system, "Observer");
    liquid_temperature_observer.generateParticles<ObserverParticles>(liquid_temperature_observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation tank_inner(tank);
    InnerRelation water_inner(water_block);
    ContactRelation water_tank_contact(water_block, {&tank});
    ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
        // SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_tank_to_vtp(tank);
        // BodyStatesRecordingToVtp write_water_to_vtp(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_tank_particle_reload_files(tank);
        // ReloadParticleIO write_water_particle_reload_files(water_block);

        /** A  Physics relaxation step. */
        // relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
        relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_tank_particles.exec(0.25);
        // random_water_particles.exec(0.25);
        tank_relaxation_step_inner.SurfaceBounding().exec();
        // water_relaxation_step_inner.SurfaceBounding().exec();
        write_tank_to_vtp.writeToFile(0);
        // write_water_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            tank_relaxation_step_inner.exec();
            // water_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
                // write_water_to_vtp.writeToFile(ite_p);
                write_tank_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of finish !" << std::endl;
        /** Output results. */
        write_tank_particle_reload_files.writeToFile(0);
        // write_water_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_inner, water_tank_contact);
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromSubShapeAndOp> tank_inner_normal_direction(tank, "InnerWall");


    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    VariableGravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity_to_water(water_block, gravity);
    
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> water_pressure_relaxation(water_inner, water_tank_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> water_density_relaxation(water_inner, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_water_density(water_inner, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force_water(water_inner, water_tank_contact);
    
    /*DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec2d, FixedDampingRate>>>
        fluid_damping(0.2, ConstructorArgs(water_inner, "Velocity", mu_water), ConstructorArgs(water_tank_contact, "Velocity", mu_water));*/
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output body states for visualization. */
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system);
    RestartIO restart_io(sph_system);

     /** WaveProbes. */
    BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("ProbeS1"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S1(probe_s1, "FreeSurfaceHeight_S1");
    BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S2(probe_s2, "FreeSurfaceHeight_S2");
    BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("ProbeS3"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S3(probe_s3, "FreeSurfaceHeight_S3");

    ReducedQuantityRecording<TotalMechanicalEnergy> write_water_mechanical_energy(water_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    tank_inner_normal_direction.exec();
    constant_gravity_to_water.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile(0);
    wave_probe_S1.writeToFile(0);
    wave_probe_S2.writeToFile(0);
    wave_probe_S3.writeToFile(0);
    write_real_body_states_to_vtp.addToWrite<Vecd>(tank, "NormalDirection"); 
    write_water_mechanical_energy.writeToFile(0);
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_tank_contact.updateConfiguration();
    }

    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 3.0; /**< End time. */
    Real output_interval = end_time / 50.0;
    Real D_Time = 0.1;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
 
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            Real Dt = get_water_advection_time_step_size.exec();
            update_water_density.exec();
            constant_gravity_to_water.exec();
            viscous_force_water.exec();
            Real relaxation_time = 0.0;
            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();

            while (relaxation_time < Dt)
            {
                dt = SMIN(get_water_time_step_size.exec(), Dt);
                //fluid_damping.exec(dt);
                /** Fluid relaxation and force computation. */
                water_pressure_relaxation.exec(dt);
                water_pressure_relaxation.exec(dt);
               
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            time_instance = TickCount::now();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList(); // water particle motion is small
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        wave_probe_S1.writeToFile();
        wave_probe_S2.writeToFile();
        wave_probe_S3.writeToFile();
        write_real_body_states_to_vtp.writeToFile();
        write_water_mechanical_energy.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}
