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

std::string probe_s1_shape = "./input/ProbeS1.stl";
std::string probe_s2_shape = "./input/ProbeS2.stl";
std::string probe_s3_shape = "./input/ProbeS3.stl";
std::string fuel_tank_outer = "./input/tank_outer_small.stl";
std::string fuel_tank_inner = "./input/tank_inner_small.stl";
std::string air = "./input/gas_small.stl";
std::string water = "./input/water_small.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real particle_spacing_ref = 0.01; /**< Reference particle spacing. */
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3, -0.3), Vec3d(0.3, 1.0, 0.3));

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 424.7;                    /*Fluid density*/
Real rho0_a = 2.668;                    /*Air density*/
Real gravity_g = 9.81;                  /*Gravity force of fluid*/
Real U_f = 2.0 * sqrt(gravity_g * 0.5); /**< Characteristic velocity. */
Real U_g = 2.0 * sqrt(gravity_g * 0.5); /**< dispersion velocity in shallow water. */
Real U_max = SMAX(U_f, U_g);
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
Real f = 1.0;
Real a = 0.01;
Real c_p_water = 3.4267e3;
Real c_p_air = 1.054e3;
Real k_water = 0.1846;
Real k_air = 0.01221;
Real diffusion_coff_water = k_water / (c_p_water * rho0_f);
Real diffusion_coff_air = k_air / (c_p_air * rho0_a);
Real mu_f = 121.79e-6;
Real mu_a = 8.73e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.0, 0);
Real initial_temperature_water = 110.15;
Real initial_temperature_air = 130.15;
std::string diffusion_species_name = "Phi";
    //----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
/*shapes from stl files*/
class TankShape : public ComplexShape
{
  public:
    explicit TankShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */

        add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale, "OuterWall");
        subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale, "InnerWall");
    }
};
class FluidShape : public ComplexShape
{
  public:
    explicit FluidShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(water, translation, length_scale);
    }
};

class AirShape : public ComplexShape
{
  public:
    explicit AirShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(air, translation, length_scale);
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
        if (time >= 0.0)
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

class ThermoWaterBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermoWaterBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>("Phi")),
          heat_flux_contact_(particles_->registerStateVariable<Real>("HeatFluxContact")),
          heat_flux_inner_(particles_->registerStateVariable<Real>("HeatFluxInner"))
    {
        this->particles_->template addEvolvingVariable<Real>("HeatFluxContact");
        this->particles_->template addVariableToWrite<Real>("HeatFluxContact");

        this->particles_->template addEvolvingVariable<Real>("HeatFluxInner");
        this->particles_->template addVariableToWrite<Real>("HeatFluxInner");
    };

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_water;
    };

  protected:
    Real *phi_;
    Real *heat_flux_inner_;
    Real *heat_flux_contact_;
};

class ThermoAirBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermoAirBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>("Phi")),
          heat_flux_contact_(particles_->registerStateVariable<Real>("HeatFluxContact")),
          heat_flux_inner_(particles_->registerStateVariable<Real>("HeatFluxInner"))
    {
        this->particles_->template addEvolvingVariable<Real>("HeatFluxContact");
        this->particles_->template addVariableToWrite<Real>("HeatFluxContact");
        this->particles_->template addEvolvingVariable<Real>("HeatFluxInner");
        this->particles_->template addVariableToWrite<Real>("HeatFluxInner");
    };

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_air;
    };

  protected:
    Real *phi_;
    Real *heat_flux_inner_;
    Real *heat_flux_contact_;
};

using HeatExchangeComplex = HeatExchangeDiffusionComplex<KernelGradientInner, KernelGradientContact, HeatIsotropicDiffusion, HeatIsotropicDiffusion>;

StdVec<Vecd> LiquidTemperatureObserverParticle()
{
    StdVec<Vecd> observation_points;

    observation_points.push_back(Vecd(0.0, 0.35, 0.0));
    observation_points.push_back(Vecd(0.0, 0.5, 0.0));
    observation_points.push_back(Vecd(0.0, 0.6, 0.0));

    return observation_points;
};

StdVec<Vecd> GasTemperatureObserverParticle()
{
    StdVec<Vecd> observation_points;

    observation_points.push_back(Vecd(0.0, 0.4, 0.0));
    observation_points.push_back(Vecd(0.0, 0.5, 0.0));
    observation_points.push_back(Vecd(0.0, 0.65, 0.0));

    return observation_points;
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
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity, HeatIsotropicDiffusion>
        (ConstructArgs(rho0_f, c_f), mu_f, ConstructArgs(diffusion_species_name, k_water, rho0_f, c_p_water));
    //water_block.generateParticles<BaseParticles, Lattice>();
    //water_block.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody air_block(sph_system, makeShared<AirShape>("AirBody"));
    air_block.defineClosure<WeaklyCompressibleFluid, Viscosity, HeatIsotropicDiffusion>
        (ConstructArgs(rho0_a, c_f), mu_a, ConstructArgs(diffusion_species_name, k_air, rho0_a, c_p_air));
    //air_block.generateParticles<BaseParticles, Lattice>();
    //air_block.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? air_block.generateParticles<BaseParticles, Reload>(air_block.getName())
        : air_block.generateParticles<BaseParticles, Lattice>();

    SolidBody tank(sph_system, makeShared<TankShape>("Tank"));
    tank.defineMaterial<Solid>();
    //tank.generateParticles<BaseParticles, Lattice>();
    //tank.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? tank.generateParticles<BaseParticles, Reload>(tank.getName())
        : tank.generateParticles<BaseParticles, Lattice>();

    ObserverBody gas_temperature_observer(sph_system, "GasTemperatureObserver");
    gas_temperature_observer.generateParticles<ObserverParticles>(GasTemperatureObserverParticle());

    ObserverBody liquid_temperature_observer(sph_system, "LiquidTemperatureObserver");
    liquid_temperature_observer.generateParticles<ObserverParticles>(LiquidTemperatureObserverParticle());

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
    ContactRelation water_tank_contact(water_block, {&tank});
    ContactRelation water_air_tank_contact(water_block, {&air_block, &tank});

    InnerRelation air_block_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_tank_contact(air_block, {&tank});
    ContactRelation air_water_tank_contact(air_block, {&water_block, &tank});

    InnerRelation tank_inner(tank);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_complex(water_block_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_complex(air_block_inner, {&air_water_contact, &air_tank_contact});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
        SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
        SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_tank_to_vtp(tank);
        BodyStatesRecordingToVtp write_air_to_vtp(air_block);
        BodyStatesRecordingToVtp write_water_to_vtp(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_tank_particle_reload_files(tank);
        ReloadParticleIO write_air_particle_reload_files(air_block);
        ReloadParticleIO write_water_particle_reload_files(water_block);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner_tank(tank_inner);
        RelaxationStepInner relaxation_step_inner_air(air_block_inner);
        RelaxationStepInner relaxation_step_inner_water(water_block_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_tank_particles.exec(0.25);
        random_air_particles.exec(0.25);
        random_water_particles.exec(0.25);
        relaxation_step_inner_tank.SurfaceBounding().exec();
        relaxation_step_inner_air.SurfaceBounding().exec();
        relaxation_step_inner_water.SurfaceBounding().exec();
        write_tank_to_vtp.writeToFile(0);
        write_air_to_vtp.writeToFile(0);
        write_water_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner_tank.exec();
            relaxation_step_inner_air.exec();
            relaxation_step_inner_water.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_tank_to_vtp.writeToFile(ite_p);
                write_air_to_vtp.writeToFile(ite_p);
                write_water_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_tank_particle_reload_files.writeToFile(0);
        write_air_particle_reload_files.writeToFile(0);
        write_water_particle_reload_files.writeToFile(0);
        return 0;
    };

    ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, {&water_block});
    ContactRelation gas_temperature_observer_contact(gas_temperature_observer, {&air_block});

    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    //GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    //CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    //SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    //SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
    //SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(tank);

    VariableGravity gravity(Vecd(0.0, -gravity_g , 0.0));
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

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_water_advection_time_step_size(water_block, U_max, 0.15);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_air_advection_time_step_size(air_block, U_max, 0.15);

    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block, 0.3);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block, 0.3);
    
    // Define diffusion coefficient
    HeatIsotropicDiffusion water_heat_diffusion("Phi", "Phi", k_water, rho0_f, c_p_water);
    HeatIsotropicDiffusion air_heat_diffusion("Phi", "Phi", k_air, rho0_a, c_p_air);

    Dynamics1Level<HeatExchangeComplex> water_heat_exchange_complex(water_block_inner, water_air_contact,  &air_heat_diffusion);
    Dynamics1Level<HeatExchangeComplex> air_heat_exchange_complex(air_block_inner, air_water_contact,  &water_heat_diffusion);

    SimpleDynamics<ThermoWaterBodyInitialCondition> water_diffusion_initial_condition(water_block);
    SimpleDynamics<ThermoAirBodyInitialCondition> air_diffusion_initial_condition(air_block);

    GetDiffusionTimeStepSize get_diffusion_time_step_size_water(water_block, &water_heat_diffusion);
    GetDiffusionTimeStepSize get_diffusion_time_step_size_air(air_block, &air_heat_diffusion);

    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting water_particle_sorting(water_block);
    ParticleSorting air_particle_sorting(air_block);
   
    /** WaveProbes. */
    BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("ProbeS1"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S1(probe_s1, "FreeSurfaceHeight_S1", 1);
    BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S2(probe_s2, "FreeSurfaceHeight_S2", 1);
    BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("ProbeS3"));
    ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
        wave_probe_S3(probe_s3, "FreeSurfaceHeight_S3", 1);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure"); 
    body_states_recording.addToWrite<Real>(air_block, "Pressure");          // output for debug
    body_states_recording.addToWrite<Vecd>(water_block, "Velocity");          // output for debug
    body_states_recording.addToWrite<Vecd>(air_block, "Velocity");            // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<Real>(air_block, "Density");             // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(tank, "NormalDirection"); // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Phi");
    body_states_recording.addToWrite<Real>(air_block, "Phi");  
    RestartIO restart_io(sph_system);
   
    ReducedQuantityRecording<TotalMechanicalEnergy> write_water_mechanical_energy(water_block, gravity);
    ReducedQuantityRecording<TotalMechanicalEnergy> write_air_mechanical_energy(air_block, gravity);

    ObservedQuantityRecording<Real> write_temperature_liquid("Phi", liquid_temperature_observer_contact);
    ObservedQuantityRecording<Real> write_temperature_gas("Phi", gas_temperature_observer_contact);

    ReducedQuantityRecording<QuantitySummation<Real>> write_water_heat_flux_inner(water_block, "HeatFluxInner");
    ReducedQuantityRecording<QuantitySummation<Real>> write_air_heat_flux_inner(air_block, "HeatFluxInner");
    ReducedQuantityRecording<QuantitySummation<Real>> write_water_heat_flux_contact(water_block, "HeatFluxContact");
    ReducedQuantityRecording<QuantitySummation<Real>> write_air_heat_flux_contact(air_block, "HeatFluxContact");
    
    ReducedQuantityRecording<Average<QuantitySummation<Real>>> water_everage_temperature(water_block, "Phi");
    ReducedQuantityRecording<Average<QuantitySummation<Real>>> air_everage_temperature(air_block, "Phi");
    ReducedQuantityRecording<QuantityMax<Real>> water_max_temperature(water_block, "Phi");
    ReducedQuantityRecording<QuantityMax<Real>> air_max_temperature(air_block, "Phi");
    
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    
    free_stream_surface_indicator.exec();
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
    water_diffusion_initial_condition.exec();
    air_diffusion_initial_condition.exec();
    wave_probe_S1.writeToFile();
    wave_probe_S2.writeToFile();
    wave_probe_S3.writeToFile();
    write_water_mechanical_energy.writeToFile(0);
    write_air_mechanical_energy.writeToFile(0);
    write_temperature_liquid.writeToFile(0);
    write_temperature_gas.writeToFile(0);
    write_water_heat_flux_inner.writeToFile(0);
    write_air_heat_flux_inner.writeToFile(0);
    write_water_heat_flux_contact.writeToFile(0);
    write_air_heat_flux_contact.writeToFile(0);
    water_everage_temperature.writeToFile(0);
    air_everage_temperature.writeToFile(0);
    water_max_temperature.writeToFile(0);
    air_max_temperature.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 10.0;
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
            Real dt_water(0.0), dt_air(0.0);

            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                Real dt_thermal_water = get_diffusion_time_step_size_water.exec();
                Real dt_thermal_air = get_diffusion_time_step_size_air.exec();
                dt = SMIN(SMIN(dt_f, dt_thermal_water), SMIN(dt_thermal_air, dt_a), Dt);
                //dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                if (physical_time >= 0.0)
                {
                    water_heat_exchange_complex.exec(dt);
                    air_heat_exchange_complex.exec(dt);
                }

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                dt_water = dt_f;
                dt_air = dt_a;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << " Dt_water = " << Dt_f << " Dt_air = " << Dt_a
                          << " dt_f = " << dt_water << " dt_a = " << dt_air<< "\n";

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
            air_block_inner.updateConfiguration();

            water_complex.updateConfiguration();
            air_complex.updateConfiguration();
            free_stream_surface_indicator.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            if (physical_time >= 0.0)
            {
                wave_probe_S1.writeToFile();
                wave_probe_S2.writeToFile();
                wave_probe_S3.writeToFile();
                write_water_mechanical_energy.writeToFile();
                write_air_mechanical_energy.writeToFile();
                write_temperature_liquid.writeToFile();
                write_temperature_gas.writeToFile();
                write_water_heat_flux_inner.writeToFile();
                write_air_heat_flux_inner.writeToFile();
                write_water_heat_flux_contact.writeToFile();
                write_air_heat_flux_contact.writeToFile();
                water_everage_temperature.writeToFile();
                air_everage_temperature.writeToFile();
                water_max_temperature.writeToFile();
                air_max_temperature.writeToFile();
            }
        }
        /*wave_probe_S1.writeToFile();
        wave_probe_S2.writeToFile();
        wave_probe_S3.writeToFile();*/
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
