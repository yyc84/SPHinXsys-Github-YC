﻿/**
 * @file	water entry and exit.cpp
 * @brief	2D water entry and exit example with surface wettability considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach, 
 *          especially when coupled with the wetting.
 * @author  Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#include "wetting_coupled_spatial_temporal_method.h"

using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real insert_cylinder_radius = 0.055;                              /**< Cylinder radius. */
Real DL = 8.0 * insert_cylinder_radius;                          /**< Water tank length. */
Real DH = 7.0 * insert_cylinder_radius;                          /**< Water tank height. */
Real LL = DL;                                                     /**< Water column length. */
Real LH = 3.0 * insert_cylinder_radius;                           /**< Water column height. */
Real particle_spacing_ref = 2.0 * insert_cylinder_radius / 40.0;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                               /**< Thickness of tank wall. */
Vec2d insert_cylinder_center(0.5 * DL, LH + 0.15);                /**< Location of the cylinder center. */
Vecd tethering_point(0.5 * DL, DH);                               /**< The tethering point. */
StdVec<Vecd> observer_location = {insert_cylinder_center};        /**< Observation point. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Fluid density. */
Real rho0_s = 0.5;                       /**< Cylinder density. */
Real gravity_g = 9.81;                   /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
Real mu_f = 8.9e-7;                      /**< Water dynamics viscosity. */
//----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
Real diffusion_coeff = 330.578 * pow(particle_spacing_ref,2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;       /**< fluid moisture. */
Real cylinder_moisture = 0.0;    /**< cylinder moisture. */
Real wall_moisture = 1.0;        /**< wall moisture. */
//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block;
    water_block.push_back(Vecd(0.0, 0.0));
    water_block.push_back(Vecd(0.0, LH));
    water_block.push_back(Vecd(LL, LH));
    water_block.push_back(Vecd(LL, 0.0));
    water_block.push_back(Vecd(0.0, 0.0));

    return water_block;
}
class WettingFluidBody : public MultiPolygonShape
{
  public:
    explicit WettingFluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
class WettingFluidBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
  public:
    WettingFluidBodyMaterial()
        : DiffusionReaction<WeaklyCompressibleFluid>({"Phi"}, SharedPtr<NoReaction>(), rho0_f, c_f, mu_f)
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
    };
};
using DiffusionFluidParticles = DiffusionReactionParticles<BaseParticles, WettingFluidBodyMaterial>;
class WettingFluidBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionFluidParticles>
{
  protected:
    size_t phi_;

  public:
    explicit WettingFluidBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionFluidParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = fluid_moisture;
    };
};
//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall;
    outer_wall.push_back(Vecd(-BW, -BW));
    outer_wall.push_back(Vecd(-BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, -BW));
    outer_wall.push_back(Vecd(-BW, -BW));

    return outer_wall;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall;
    inner_wall.push_back(Vecd(0.0, 0.0));
    inner_wall.push_back(Vecd(0.0, DH));
    inner_wall.push_back(Vecd(DL, DH));
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));

    return inner_wall;
}
class WettingWallBody : public MultiPolygonShape
{
  public:
    explicit WettingWallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
class WettingWallBodyMaterial : public DiffusionReaction<Solid>
{
  public:
    WettingWallBodyMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
    };
};
using DiffusionWallParticles = DiffusionReactionParticles<SolidParticles, WettingWallBodyMaterial>;
class WettingWallBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionWallParticles>
{
  protected:
    size_t phi_;

  public:
    explicit WettingWallBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionWallParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = wall_moisture;
    };
};
//----------------------------------------------------------------------
//	Definition for cylinder body
//----------------------------------------------------------------------
class WettingCylinderBody : public MultiPolygonShape
{
  public:
    explicit WettingCylinderBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(insert_cylinder_center, insert_cylinder_radius, 100, ShapeBooleanOps::add);
    }
};
class WettingCylinderBodyMaterial : public DiffusionReaction<Solid>
{
  public:
    WettingCylinderBodyMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>(), rho0_s)
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi",diffusion_coeff);
    };
};
using DiffusionCylinderParticles = DiffusionReactionParticles<SolidParticles, WettingCylinderBodyMaterial>;
class WettingCylinderBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionCylinderParticles>
{
  protected:
    size_t phi_;

  public:
    explicit WettingCylinderBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionCylinderParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = cylinder_moisture;
    };
};

//----------------------------------------------------------------------
//	Set topology for wetting bodies
//----------------------------------------------------------------------
using CylinderFluidDiffusionDirichlet = DiffusionRelaxationDirichlet<DiffusionCylinderParticles, DiffusionFluidParticles>;
class ThermalRelaxationComplex
    : public DiffusionRelaxationRK2<
          ComplexInteraction<CylinderFluidDiffusionDirichlet>>
{
  public:
    explicit ThermalRelaxationComplex(BaseContactRelation &body_contact_relation_Dirichlet)
        : DiffusionRelaxationRK2<ComplexInteraction<CylinderFluidDiffusionDirichlet>>(body_contact_relation_Dirichlet){};
    virtual ~ThermalRelaxationComplex(){};
};
//------------------------------------------------------------------------------
// Constrained part for Simbody
//------------------------------------------------------------------------------
MultiPolygon createSimbodyConstrainShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(insert_cylinder_center, insert_cylinder_radius, 100, ShapeBooleanOps::add);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    GlobalStaticVariables::physical_time_=0.0;
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WettingFluidBody>("WaterBody"));
    water_block.defineParticlesAndMaterial<DiffusionFluidParticles, WettingFluidBodyMaterial>();
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WettingWallBody>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<DiffusionWallParticles, WettingWallBodyMaterial>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

    SolidBody cylinder(sph_system, makeShared<WettingCylinderBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 1.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineParticlesAndMaterial<DiffusionCylinderParticles, WettingCylinderBodyMaterial>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
        : cylinder.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticleGenerator>(observer_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ComplexRelation water_block_complex(water_block_inner, {&wall_boundary, &cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    ContactRelation fluid_observer_contact(fluid_observer, {&cylinder});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, {&cylinder});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&cylinder});
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block, gravity_ptr);
    InteractionWithUpdate<fluid_dynamics::WettingCoupledSpatialTemporalFreeSurfaceIdentificationComplex>
        free_stream_surface_indicator(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> fluid_density_by_summation(water_block_complex);
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("SurfaceIndicator");
    cylinder.addBodyStateForRecording<Real>("Density");
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> fluid_pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> fluid_density_relaxation(water_block_complex);   
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction(water_block_complex);
    //----------------------------------------------------------------------
    //	Define the wetting diffusion dynamics used in the simulation.
    //----------------------------------------------------------------------
    SimpleDynamics<WettingFluidBodyInitialCondition> Wetting_water_initial_condition(water_block);
    SimpleDynamics<WettingWallBodyInitialCondition> Wetting_wall_initial_condition(wall_boundary);
    SimpleDynamics<WettingCylinderBodyInitialCondition> Wetting_cylinder_initial_condition(cylinder);
    GetDiffusionTimeStepSize<DiffusionCylinderParticles> get_thermal_time_step(cylinder);
    ThermalRelaxationComplex thermal_relaxation_complex(cylinder_contact);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> fluid_viscous_force_on_inserted_body(cylinder_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> 
        fluid_pressure_force_on_inserted_body(cylinder_contact, fluid_viscous_force_on_inserted_body);    
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass proeprties of the fixed spot. */
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    SolidBodyPartForSimbody cylinder_constraint_area(cylinder, makeShared<MultiPolygonShape>(createSimbodyConstrainShape(cylinder), "cylinder"));
    /** Mass properties of the consrained spot. */
    SimTK::Body::Rigid tethered_spot_info(*cylinder_constraint_area.body_part_mass_properties_);
    /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], 0.0)),
                                          fixed_spot_info, SimTK::Transform(SimTKVec3(0)));
    /** Mobility of the tethered spot.
     * Set the mass center as the origin location of the planar mobilizer
     */
    Vecd displacement0 = cylinder_constraint_area.initial_mass_center_ - tethering_point;
    SimTK::MobilizedBody::Planar tethered_spot(fixed_spot,
                                               SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], 0.0)), tethered_spot_info, SimTK::Transform(SimTKVec3(0)));
    // discreted forces acting on the bodies
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-9.81), 0.0), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));
    SimTK::State state = MBsystem.realizeTopology();

    /** Time steping method for multibody system.*/
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH..
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RestartIO restart_io(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_cylinder_displacement("Position", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    cylinder_normal_direction.exec();
    Wetting_water_initial_condition.exec();
    Wetting_wall_initial_condition.exec();
    Wetting_cylinder_initial_condition.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    free_stream_surface_indicator.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_block_complex.updateConfiguration();
        fluid_observer_contact.updateConfiguration();
        cylinder.updateCellLinkedList();
        water_block_inner.updateConfiguration();
        cylinder_inner.updateConfiguration();
        cylinder_contact.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 0.7;
    Real output_interval = end_time/70.0;
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
    write_cylinder_displacement.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_step_initialization.exec();
            Real Dt = fluid_advection_time_step.exec();
            
            fluid_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            fluid_viscous_force_on_inserted_body.exec();
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                fluid_pressure_relaxation.exec(dt);
                fluid_pressure_force_on_inserted_body.exec();
                fluid_density_relaxation.exec(dt);
                thermal_relaxation_complex.exec(dt);

                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                constraint_tethered_spot.exec();
               
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_cylinder_displacement.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedListWithParticleSort(100);
            cylinder.updateCellLinkedList();
            water_block_inner.updateConfiguration();
            cylinder_inner.updateConfiguration();
            cylinder_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
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

    if (sph_system.generate_regression_data_)
    {
        write_cylinder_displacement.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_cylinder_displacement.testResult();
    }

    return 0;
};
