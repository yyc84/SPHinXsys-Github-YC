/**
 * @file	water entry and exit.cpp
 * @brief	3D water entry and exit example with surface wetting considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach,
 *          especially when coupled with the wetting.
 * @author  Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real particle_spacing_ref = 0.005; /**< Reference particle spacing. */
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3, -0.3), Vec3d(0.3, 0.5, 0.3));
Real Tank_L = 0.2;                          /**< Water tank length. */
Real Tank_H = 0.35;                          /**< Water tank height. */
Real Tank_W = 0.2;                          /**< Water tank width. */
Real Liquid_L = 0.2;                                             /**< Water column length. */
Real Liquid_H = 0.175;                          /**< Water column height. */
Real Liquid_W = 0.2;                                             /**< Water column length. */
Real Gas_L = 0.2;                            /**< Water column length. */
Real Gas_H = 0.175;                          /**< Water column height. */
Real Gas_W = 0.2;
Real BW = particle_spacing_ref * 4;                       /**< Thickness of tank wall. */

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                       /**< Fluid density. */
Real gravity_g = 9.81;                   /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * Liquid_L); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
Real mu_f = 653.9e-6;                    /**< Water dynamics viscosity. */
Real rho0_a = 1.226; /**< Reference density of air. */
Real mu_a = 20.88e-6;


Real f = 1.63;
Real a = 0.01;

//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
class FluidShape : public ComplexShape
{
  public:
    explicit FluidShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * Liquid_L, 0.5 * Liquid_H, 0.5 * Liquid_W);
        Transform translation_water(Vecd(0.0, -0.5 * Liquid_H, 0.0));
        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
    }
};

//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
class TankShape : public ComplexShape
{
  public:
    explicit TankShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * Tank_L + BW, 0.5 * Tank_H + BW, 0.5 * Tank_W + BW);
        Vecd halfsize_inner(0.5 * Tank_L, 0.5 * Tank_H, 0.5 * Tank_W);
        Transform translation_wall(Vecd(0.0, 0.0, 0.0));
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
    }
};

class AirShape : public ComplexShape
{
  public:
    explicit AirShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_inner(0.5 * Gas_L, 0.5 * Gas_H, 0.5 * Gas_W);
        Transform translation_wall(Vecd(0.0, 0.5 * Gas_H, 0.0));
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);

    }
};


//----------------------------------------------------------------------
//	The diffusion model of wetting
//----------------------------------------------------------------------
using CylinderFluidDiffusionDirichlet =
    DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;

class VariableGravity : public Gravity
{

  public:
    VariableGravity(Vecd gravity_vector) : Gravity(gravity_vector) {};
    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        Real time = physical_time;
        Vecd acceleration = reference_acceleration_;
        if (time >= 1.0)
        {
            acceleration[0] = -4.0 * PI * PI * f * f * a * sin(2 * PI * f * (time - 1));
        }
        // global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
        return acceleration;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    //BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -BW), Vec3d(DL + BW, DW + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<FluidShape>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody air_block(sph_system, makeShared<AirShape>("AirBody"));
    air_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_a, c_f), mu_a);
    air_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<TankShape>("Tank"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

   /* SolidBody cylinder(sph_system, makeShared<WettingCylinderBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 2.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();*/

   /* ObserverBody cylinder_observer(sph_system, "CylinderObserver");
    cylinder_observer.generateParticles<ObserverParticles>(observer_location);

    ObserverBody wetting_observer(sph_system, "WettingObserver");
    wetting_observer.generateParticles<ObserverParticles>(wetting_observer_location);*/
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation water_air_wall_contact(water_block, {&air_block, &wall_boundary});

    InnerRelation air_block_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_wall_contact(air_block, {&wall_boundary});
    ContactRelation air_water_wall_contact(air_block, {&water_block, &wall_boundary});

    /*InnerRelation cylinder_inner(cylinder);
    ContactRelation cylinder_contact_water_only(cylinder, {&water_block});
    ContactRelation cylinder_contact(cylinder, {&water_block, &air_block});*/

    /*ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder});
    ContactRelation wetting_observer_contact(wetting_observer, {&cylinder});*/
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_complex(water_block_inner, {&water_air_contact, &water_wall_contact});
    ComplexRelation air_complex(air_block_inner, {&air_water_contact, &air_wall_contact});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    //if (sph_system.RunParticleRelaxation())
    //{
    //    /** body topology only for particle relaxation */
    //    InnerRelation cylinder_inner(cylinder);
    //    //----------------------------------------------------------------------
    //    //	Methods used for particle relaxation.
    //    //----------------------------------------------------------------------
    //    using namespace relax_dynamics;
    //    SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
    //    /** Write the body state to Vtp file. */
    //    BodyStatesRecordingToVtp write_inserted_body_to_vtp(cylinder);
    //    /** Write the particle reload files. */
    //    ReloadParticleIO write_particle_reload_files(cylinder);
    //    /** A  Physics relaxation step. */
    //    RelaxationStepInner relaxation_step_inner(cylinder_inner);
    //    //----------------------------------------------------------------------
    //    //	Particle relaxation starts here.
    //    //----------------------------------------------------------------------
    //    random_inserted_body_particles.exec(0.25);
    //    relaxation_step_inner.SurfaceBounding().exec();
    //    write_inserted_body_to_vtp.writeToFile(0);
    //    //----------------------------------------------------------------------
    //    //	Relax particles of the insert body.
    //    //----------------------------------------------------------------------
    //    int ite_p = 0;
    //    while (ite_p < 1000)
    //    {
    //        relaxation_step_inner.exec();
    //        ite_p += 1;
    //        if (ite_p % 200 == 0)
    //        {
    //            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    //            write_inserted_body_to_vtp.writeToFile(ite_p);
    //        }
    //    }
    //    std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    //    /** Output results. */
    //    write_particle_reload_files.writeToFile(0);
    //    return 0;
    //}
    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    //GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    //CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    //SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    //SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
    //SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    //SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    VariableGravity gravity(Vecd(0.0, -gravity_g , 0.0));
    SimpleDynamics<GravityForce<VariableGravity>> constant_gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce<VariableGravity>> constant_gravity_to_air(air_block, gravity);
    InteractionDynamics<fluid_dynamics::BoundingFromWall> air_near_wall_bounding(air_wall_contact);
    //InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        water_pressure_relaxation(water_block_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        water_density_relaxation(water_block_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        air_pressure_relaxation(air_block_inner, air_water_contact, air_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        air_density_relaxation(air_block_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface>
        update_water_density_by_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_block_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_block_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<BulkParticles>>
        water_transport_correction(water_block_inner, water_air_contact, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall>
        water_viscous_acceleration(water_block_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall>
        air_viscous_acceleration(air_block_inner, air_water_contact, air_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_water_advection_time_step_size(water_block, U_max, 0.15);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_air_advection_time_step_size(air_block, U_max, 0.15);

    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block, 0.3);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block, 0.3);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    //InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    //InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(water_density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting water_particle_sorting(water_block);
    ParticleSorting air_particle_sorting(air_block);
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    //SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    //SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    //SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the fixed spot. */
    //GeometricShapeBall fix_spot_shape(cylinder_center, cylinder_radius);
    //StructureSystemForSimbody ball_multibody(cylinder, fix_spot_shape);
    //SimTK::Body::Rigid ball_info(*ball_multibody.body_part_mass_properties_);
    /** Mobility of the tethered spot.
     * Set the mass center as the origin location of the planar mobilizer
     */
    //SimTK::MobilizedBody::Free ball_mob(matter.Ground(), SimTK::Transform(SimTKVec3(G[0], G[1], G[2])),ball_info, SimTK::Transform(SimTKVec3(0)));

    // discrete forces acting on the bodies
    //SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, 0.0, Real(-gravity_g)), 0.0);
    //SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    //SimTK::State state = MBsystem.realizeTopology();

    /** Time stepping method for multibody system.*/
    //SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    //integ.setAccuracy(1e-3);
    //integ.setAllowInterpolation(false);
    //integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH..
    //----------------------------------------------------------------------
    //ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody> force_on_ball(ball_multibody, MBsystem, ball_mob, integ);
    //SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody> constraint_tethered_spot(ball_multibody, MBsystem, ball_mob, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");          // output for debug
   // body_states_recording.addToWrite<Vecd>(cylinder, "Velocity");          // output for debug
    body_states_recording.addToWrite<Vecd>(water_block, "Velocity");          // output for debug
   // body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection"); // output for debug
    RestartIO restart_io(sph_system);
    //ObservedQuantityRecording<Vecd> write_cylinder_displacement("Position", cylinder_observer_contact);
    //ObservedQuantityRecording<Real> write_cylinder_wetting("Phi", wetting_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //cylinder_normal_direction.exec();
    //wetting_water_initial_condition.exec();
   // wetting_wall_initial_condition.exec();
    //wetting_cylinder_initial_condition.exec();
    //Real dt_thermal = get_thermal_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 11.0;
    Real output_interval = end_time / 70.0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //write_cylinder_displacement.writeToFile(number_of_iterations);
    //write_cylinder_wetting.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            constant_gravity_to_water.exec();
            constant_gravity_to_air.exec();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_f, Dt_a);

            update_water_density_by_summation.exec();
            update_air_density_by_summation.exec();
            water_viscous_acceleration.exec();
            air_viscous_acceleration.exec();
            air_transport_correction.exec();
            air_near_wall_bounding.exec();
            water_transport_correction.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            //viscous_force_from_fluid.exec();
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                //pressure_force_from_fluid.exec();

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                //cylinder_wetting.exec(dt);

                //integ.stepBy(dt);
                //SimTK::State &state_for_update = integ.updAdvancedState();
                //force_on_bodies.clearAllBodyForces(state_for_update);
                //force_on_bodies.setOneBodyForce(state_for_update, ball_mob, force_on_ball.exec());
                //constraint_tethered_spot.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                /*if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_cylinder_displacement.writeToFile(number_of_iterations);
                    write_cylinder_wetting.writeToFile(number_of_iterations);
                }*/
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                water_particle_sorting.exec();
                air_particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            air_block.updateCellLinkedList();
            //cylinder.updateCellLinkedList();

            water_block_inner.updateConfiguration();
            air_block_inner.updateConfiguration();
            //cylinder_inner.updateConfiguration();
            //cylinder_contact.updateConfiguration();

            water_complex.updateConfiguration();
            air_complex.updateConfiguration();
            free_stream_surface_indicator.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    /* if (sph_system.GenerateRegressionData())
    {
        write_cylinder_displacement.generateDataBase(1.0e-3);
        write_cylinder_wetting.generateDataBase(1.0e-3);
    }
    else
    {
        write_cylinder_displacement.testResult();
        write_cylinder_wetting.testResult();
    }*/

    return 0;
};
