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
Real cylinder_radius = 0.3;                             /**< Cylinder radius. */
Real DL = 8.0 * cylinder_radius;                          /**< Water tank length. */
Real DH = 7.0 * cylinder_radius;                          /**< Water tank height. */
Real DW = 8.0 * cylinder_radius;                          /**< Water tank width. */
Real LL = DL;                                             /**< Water column length. */
Real LH = 3.0 * cylinder_radius;                          /**< Water column height. */
Real LW = DW;                                             /**< Water column length. */
Real particle_spacing_ref = 2.0 * cylinder_radius / 20.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                       /**< Thickness of tank wall. */
Vec3d cylinder_center(0.5 * DL, 0.5 * DW, DH - 2.0*cylinder_radius); /**< Location of the cylinder center. */
Vecd tethering_point(0.5 * DL, 0.5 * DW, DH);             /**< The tethering point. */
StdVec<Vecd> observer_location = {cylinder_center};       /**< Displacement observation point. */
StdVec<Vecd> wetting_observer_location ={cylinder_center - Vecd(0.0, 0.0, cylinder_radius)}; /**< wetting observation point. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
//Real rho0_f = 1000.0;                       /**< Fluid density. */
//Real rho0_s = 700;    /**< Cylinder density. */
//Real BallVolume = 4.0 / 3.0 * Pi * cylinder_radius * cylinder_radius * cylinder_radius;
//Real BallMass = rho0_s * BallVolume;
//Real gravity_g = 9.81;                   /**< Gravity. */
//Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
//Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
//Real mu_f = 653.9e-6;                    /**< Water dynamics viscosity. */
//Real rho0_a = 1.226; /**< Reference density of air. */
//Real mu_a = 20.88e-6;

Real resolution_ref = 0.005; /* Initial particle spacing*/
/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3, -0.3), Vec3d(0.3, 0.5, 0.3));
/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;                     /*Fluid density*/
Real rho0_a = 1.226;                      /*Air density*/
Real gravity_g = 9.81;                    /*Gravity force of fluid*/
Real U_f = 2.0 * sqrt(gravity_g * 0.174); /**< Characteristic velocity. */
Real U_g = 2.0 * sqrt(gravity_g * 0.174); /**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);         /**< Reference sound speed. */
Real f = 1.63;
Real a = 0.0075;
Real mu_water = 653.9e-6;
Real mu_air = 20.88e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0, 0);

//std::string fuel_tank_outer = "./input/validation_tank_outer_slim.stl";
//std::string fuel_tank_inner = "./input/validation_tank_inner.stl";
//std::string water_05 = "./input/validation_water.stl";
//std::string air_05 = "./input/validation_air.stl";
std::string probe_s1_shape = "./input/ProbeS1.stl";
std::string probe_s2_shape = "./input/ProbeS2.stl";
std::string probe_s3_shape = "./input/ProbeS3.stl";

std::string fuel_tank_outer = "./input/square_tank_outer.stl";
std::string fuel_tank_inner = "./input/square_tank_inner.stl";
std::string water_05 = "./input/square_tank_water.stl";
std::string air_05 = "./input/square_tank_gas.stl";

    //----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 100.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real cylinder_moisture = 0.0;                                /**< cylinder moisture. */
Real wall_moisture = 1.0;                                    /**< wall moisture. */
//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
//class WettingFluidBody : public ComplexShape
//{
//  public:
//    explicit WettingFluidBody(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        Vecd halfsize_water(0.5 * LL, 0.5 * LW, 0.5 * LH);
//        Transform translation_water(halfsize_water);
//        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
//    }
//};

//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
//class WettingWallBody : public ComplexShape
//{
//  public:
//    explicit WettingWallBody(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DW + BW, 0.5 * DH + BW);
//        Vecd halfsize_inner(0.5 * DL, 0.5 * DW, 0.5 * DH);
//        Transform translation_wall(halfsize_inner);
//        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
//        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
//    }
//};

//class AirBody : public ComplexShape
//{
//  public:
//    explicit AirBody(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        Vecd halfsize_inner(0.5 * DL, 0.5 * DW, 0.5 * DH);
//        Transform translation_wall(halfsize_inner);
//        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
//
//        Vecd halfsize_water(0.5 * LL, 0.5 * LW, 0.5 * LH);
//        Transform translation_water(halfsize_water);
//        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
//
//        subtract<GeometricShapeBall>(cylinder_center, cylinder_radius);
//    }
//};

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

  public:
    VariableGravity(Vecd gravity_vector) : Gravity(gravity_vector) {};
    Vecd InducedAcceleration(const Vecd &position, Real physical_time) const
    {
        Real time = physical_time;
        Vecd acceleration = reference_acceleration_;
        if (time >= 2.0)
        {
            acceleration[0] = -4.0 * PI * PI * f * f * a * sin(2 * PI * f * (time - 2));
        }
        // global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
        return acceleration;
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

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    //BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -BW), Vec3d(DL + BW, DW + BW, DH + BW));
    //SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    system.setRunParticleRelaxation(true);
    system.setReloadParticles(false);
    system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_water);
    water_block.defineBodyLevelSetShape()->writeLevelSet(system);
    //water_block.generateParticles<BaseParticles, Lattice>();
    /*(!system.RunParticleRelaxation() && system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();*/

    FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
    air_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_a, c_f), mu_air);
    air_block.generateParticles<BaseParticles, Lattice>();
    //air_block.defineBodyLevelSetShape()->writeLevelSet(system);
    /*(!system.RunParticleRelaxation() && system.ReloadParticles())
        ? air_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : air_block.generateParticles<BaseParticles, Lattice>();*/

    SolidBody tank(system, makeShared<Tank>("Tank"));
    tank.defineMaterial<Solid>();
    tank.defineBodyLevelSetShape()->writeLevelSet(system);
   // tank.generateParticles<BaseParticles, Lattice>();
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? tank.generateParticles<BaseParticles, Reload>(tank.getName())
        : tank.generateParticles<BaseParticles, Lattice>();

    

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation tank_inner(tank);
    ContactRelation tank_water_contact(tank, {&water_block});
    ContactRelation tank_air_contact(tank, {&air_block});
    ComplexRelation tank_complex(tank_inner, {&tank_water_contact, &tank_air_contact});

    InnerRelation water_block_inner(water_block);
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation water_tank_contact(water_block, {&tank});
    ContactRelation water_air_tank_contact(water_block, {&air_block, &tank});

    InnerRelation air_block_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_tank_contact(air_block, {&tank});
    ContactRelation air_water_tank_contact(air_block, {&water_block, &tank});

    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_air_tank_complex(water_block_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_water_tank_complex(air_block_inner, {&air_water_contact, &air_tank_contact});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
        //SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
        //SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_tank_to_vtp(tank);
        //BodyStatesRecordingToVtp write_water_to_vtp(water_block);
        //BodyStatesRecordingToVtp write_air_to_vtp(air_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_tank_particle_reload_files(tank);
        //ReloadParticleIO write_water_particle_reload_files(water_block);
        //ReloadParticleIO write_air_particle_reload_files(air_block);
        /** A  Physics relaxation step. */
        RelaxationStepInner tank_relaxation_step_inner(tank_inner);
        //RelaxationStepInner water_relaxation_step_inner(water_block_inner);
        //RelaxationStepInner air_relaxation_step_inner(air_block_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_tank_particles.exec(0.25);
        //random_water_particles.exec(0.25);
        //random_air_particles.exec(0.25);
        tank_relaxation_step_inner.SurfaceBounding().exec();
        //water_relaxation_step_inner.SurfaceBounding().exec();
        //air_relaxation_step_inner.SurfaceBounding().exec();
        write_tank_to_vtp.writeToFile(0);
        //write_water_to_vtp.writeToFile(0);
        //write_air_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            tank_relaxation_step_inner.exec();
            //water_relaxation_step_inner.exec();
            //air_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_tank_to_vtp.writeToFile(ite_p);
                //write_water_to_vtp.writeToFile(ite_p);
                //write_air_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_tank_particle_reload_files.writeToFile(0);
        //write_water_particle_reload_files.writeToFile(0);
        //write_air_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    //GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    //CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    //SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    //SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
   // SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

    //SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(tank);
    SimpleDynamics<NormalDirectionFromSubShapeAndOp> inner_normal_direction(tank, "InnerWall");
    //SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    VariableGravity gravity(Vecd(0.0, -gravity_g, 0.0 ));
    SimpleDynamics<GravityForce<VariableGravity>> constant_gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce<VariableGravity>> constant_gravity_to_air(air_block, gravity);
    InteractionDynamics<fluid_dynamics::BoundingFromWall> air_near_wall_bounding(air_tank_contact);
    //InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_tank_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        water_pressure_relaxation(water_block_inner, water_air_contact, water_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        water_density_relaxation(water_block_inner, water_air_contact, water_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        air_pressure_relaxation(air_block_inner, air_water_contact, air_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        air_density_relaxation(air_block_inner, air_water_contact, air_tank_contact);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface>
        update_water_density_by_summation(water_block_inner, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_block_inner, air_water_contact, air_tank_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_block_inner, air_water_contact, air_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<BulkParticles>>
        water_transport_correction(water_block_inner, water_air_contact, water_tank_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall>
        water_viscous_acceleration(water_block_inner, water_air_contact, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall>
        air_viscous_acceleration(air_block_inner, air_water_contact, air_tank_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_water_advection_time_step_size(water_block, U_f, 0.25);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_air_advection_time_step_size(air_block, U_g, 0.25);

    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block, 0.6);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block, 0.6);
    
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting water_particle_sorting(water_block);
    ParticleSorting air_particle_sorting(air_block);
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    //SimTK::MultibodySystem MBsystem;
    ///** The bodies or matter of the MBsystem. */
    //SimTK::SimbodyMatterSubsystem matter(MBsystem);
    ///** The forces of the MBsystem.*/
    //SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the fixed spot. */
   /* GeometricShapeBall fix_spot_shape(cylinder_center, cylinder_radius);
    StructureSystemForSimbody ball_multibody(cylinder, fix_spot_shape);
    SimTK::Body::Rigid ball_info(*ball_multibody.body_part_mass_properties_);*/
    /** Mobility of the tethered spot.
     * Set the mass center as the origin location of the planar mobilizer
     */
    /*SimTK::MobilizedBody::Free ball_mob(matter.Ground(),
                                               SimTK::Transform(SimTKVec3(G[0], G[1], G[2])),
                                               ball_info, SimTK::Transform(SimTKVec3(0)));*/

    // discrete forces acting on the bodies
   /* SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, 0.0, Real(-gravity_g)), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    SimTK::State state = MBsystem.realizeTopology();*/

    /** Time stepping method for multibody system.*/
    /*SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);*/
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH..
    //----------------------------------------------------------------------
    /*ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_ball(ball_multibody, MBsystem, ball_mob, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_tethered_spot(ball_multibody, MBsystem, ball_mob, integ);*/

    BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("ProbeS1"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S1(probe_s1, "FreeSurfaceHeight_S1", 1);
    BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S2(probe_s2, "FreeSurfaceHeight_S2", 1);
    BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("ProbeS3"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S3(probe_s3, "FreeSurfaceHeight_S3", 1);

    ReducedQuantityRecording<TotalMechanicalEnergy> write_water_mechanical_energy(water_block, gravity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");          // output for debug
   // body_states_recording.addToWrite<Vecd>(air_block, "Density");          // output for debug
    body_states_recording.addToWrite<Vecd>(water_block, "Velocity");          // output for debug
    //body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(tank, "NormalDirection"); // output for debug

    //ObservedQuantityRecording<Vecd> write_cylinder_displacement("Position", cylinder_observer_contact);
   // ObservedQuantityRecording<Real> write_cylinder_wetting("Phi", wetting_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    inner_normal_direction.exec();
    //cylinder_normal_direction.exec();
    //wetting_water_initial_condition.exec();
    //wetting_wall_initial_condition.exec();
    //wetting_cylinder_initial_condition.exec();
   // Real dt_thermal = get_thermal_time_step.exec();
    //free_stream_surface_indicator.exec();
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 12.0;
    Real output_interval = 0.1;
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
   // write_cylinder_displacement.writeToFile(number_of_iterations);
   // write_cylinder_wetting.writeToFile(number_of_iterations);
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

               
               // force_on_bodies.clearAllBodyForces(state_for_update);
               // force_on_bodies.setOneBodyForce(state_for_update, ball_mob, force_on_ball.exec());
               // constraint_tethered_spot.exec();

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

               /* if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_cylinder_displacement.writeToFile(number_of_iterations);
                    write_cylinder_wetting.writeToFile(number_of_iterations);
                }*/
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

            water_block_inner.updateConfiguration();
            water_air_contact.updateConfiguration();
            water_air_tank_complex.updateConfiguration();
            air_block_inner.updateConfiguration();
            air_water_contact.updateConfiguration();
            air_water_tank_complex.updateConfiguration();

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

    return 0;
};
