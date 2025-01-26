/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;
#define PI (3.14159265358979323846)

// general parameters for geometry
//Real resolution_ref = 0.05;   // particle spacing
//Real BW = resolution_ref * 4; // boundary width
Real DL = 2.0;              // tank length
Real DH = 4.0;                // tank height
Real DW = 2.0;                // tank width
Real LL = 2.0;                // liquid length
Real LH = 2.0;                // liquid height
Real LW = 2.0;                // liquid width

// for material properties of the fluid
//Real rho0_f = 1.0;
//Real gravity_g = 1.0;
//Real U_f = 2.0 * sqrt(gravity_g * LH);
//Real c_f = 10.0 * U_f;

/*tank parameters*/
Real resolution_ref = 0.008; /* Initial particle spacing*/
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

BoundingBox system_domain_bounds(Vec3d(-0.2, -0.05, -0.2), Vec3d(0.2, 1.0, 0.2));

//	define the water block shape
//class WaterBlock : public ComplexShape
//{
//  public:
//    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
//        Transform translation_water(halfsize_water);
//        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
//    }
//};
////	define the static solid wall boundary shape
//class Tank : public ComplexShape
//{
//  public:
//    explicit Tank(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        Vecd halfsize_outer(0.5 * LL + BW, 0.5 * DH + BW, 0.5 * LW + BW);
//        Vecd halfsize_inner(0.5 * LL, 0.5 * DH, 0.5 * LW);
//        Transform translation_wall(halfsize_inner);
//        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
//        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
//    }
//};

//StdVec<Vecd> createObservationPoints()
//{
//    StdVec<Vecd> observation_points;
//    observation_points.push_back(Vecd(LL, 0.01, 0.5 * LW));
//    observation_points.push_back(Vecd(LL, 0.1, 0.5 * LW));
//    observation_points.push_back(Vecd(LL, 0.2, 0.5 * LW));
//    observation_points.push_back(Vecd(LL, 0.24, 0.5 * LW));
//    observation_points.push_back(Vecd(LL, 0.252, 0.5 * LW));
//    observation_points.push_back(Vecd(LL, 0.266, 0.5 * LW));
//    return observation_points;
//};

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

class AirBlock : public ComplexShape
{
  public:
    explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
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
        if (time_ <= 5.0)
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
StdVec<Vecd> gas_temperature_observation_location = {Vecd(0.0, 0.4, 0.0), Vecd(0.0, 0.5, 0.0), Vecd(0.0, 0.65, 0.0)};
// the main program with commandline options

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    //BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(true);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(false);
    /* Tag for computation from restart files. 0: start with initial condition. */
    sph_system.setRestartStep(0);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    /*WaterBlock initial_water_block("WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();*/

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
    air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_a, c_f, mu_air);
    // air_block.defineBodyLevelSetShape()->writeLevelSet(system);
    air_block.generateParticles<BaseParticles, Lattice>();

    /*SolidBody tank(sph_system, makeShared<Tank>("Tank"));
    tank.defineMaterial<Solid>();
    tank.generateParticles<BaseParticles, Lattice>();*/

    SolidBody tank(sph_system, makeShared<Tank>("Tank"));
    //tank.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    tank.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
         ? tank.generateParticles<BaseParticles, Reload>(tank.getName())
         : tank.generateParticles<BaseParticles, Lattice>();

    ObserverBody water_observer(sph_system, "WaterObserver");
    water_observer.generateParticles<ObserverParticles>(liquid_temperature_observation_location);

    ObserverBody gas_observer(sph_system, "GasObserver");
    gas_observer.generateParticles<ObserverParticles>(gas_temperature_observation_location);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    InnerRelation tank_inner(tank);
    InnerRelation air_inner(air_block);
    
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
    ContactRelation water_tank_contact(water_block, {&tank});
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_tank_contact(air_block, {&tank});
    ContactRelation water_observer_contact(water_observer, {&water_block});
    ContactRelation air_observer_contact(gas_observer, { &air_block });
    ComplexRelation water_air_tank_complex(water_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_water_tank_complex(air_inner, {&air_water_contact, &air_tank_contact});
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    //Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    VariableGravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce> constant_gravity_to_air(air_block, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> tank_normal_direction(tank);
    InteractionDynamics<fluid_dynamics::BoundingFromWall> air_near_wall_bounding(air_tank_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> water_pressure_relaxation(water_inner, water_air_contact, water_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> water_density_relaxation(water_inner, water_air_contact, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> water_update_density_by_summation(water_inner, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_force_water(water_inner, water_air_contact, water_tank_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> air_pressure_relaxation(air_inner, air_water_contact, air_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> air_density_relaxation(air_inner, air_water_contact, air_tank_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>> air_update_density_by_summation(air_inner, air_water_contact, air_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_force_air(air_inner, air_water_contact, air_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>> air_transport_correction(air_inner, air_water_contact, air_tank_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_g);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vec3d>(tank, "NormalDirection");

    RestartIO restart_io(sph_system);

    BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("ProbeS1"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_S1(probe_s1, "FreeSurfaceHeight_S1");
    BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_S2(probe_s2, "FreeSurfaceHeight_S2");
    BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("ProbeS3"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>> wave_probe_S3(probe_s3, "FreeSurfaceHeight_S3");

    ReducedQuantityRecording<TotalMechanicalEnergy> write_water_mechanical_energy(water_block, gravity);
    ReducedQuantityRecording<TotalMechanicalEnergy> write_air_mechanical_energy(air_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    tank_normal_direction.exec();
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 10.0;
    Real output_interval = end_time / 20.0;
    Real D_Time = 0.1; /**< time stamps for output. */
    Real dt = 0.0; // default acoustic time step sizes
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    wave_probe_S1.writeToFile(0);
    wave_probe_S2.writeToFile(0);
    wave_probe_S3.writeToFile(0);
    //write_water_mechanical_energy.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            time_instance = TickCount::now();
            Real Dt_w = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_w, Dt_a);
            constant_gravity_to_water.exec();
            constant_gravity_to_air.exec();

            water_update_density_by_summation.exec();
            air_update_density_by_summation.exec();
            air_transport_correction.exec();
            //viscous_force_air.exec();
            //viscous_force_water.exec();
            air_near_wall_bounding.exec();
            
            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);
          
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
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

            water_block.updateCellLinkedListWithParticleSort(100);
            water_air_tank_complex.updateConfiguration();
            water_observer_contact.updateConfiguration();

            air_block.updateCellLinkedListWithParticleSort(100);
            air_water_tank_complex.updateConfiguration();
            air_observer_contact.updateConfiguration();
            
            write_real_body_states.writeToFile();

            interval_updating_configuration += TickCount::now() - time_instance;
        }

        //write_water_mechanical_energy.writeToFile(number_of_iterations);

        TickCount t2 = TickCount::now();
        wave_probe_S1.writeToFile();
        wave_probe_S2.writeToFile();
        wave_probe_S3.writeToFile();
        write_real_body_states.writeToFile();
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
