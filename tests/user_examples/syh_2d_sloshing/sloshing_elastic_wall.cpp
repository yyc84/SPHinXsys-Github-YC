/**
 * @file 	dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)
using namespace SPH;   // Namespace cite here.

Real particle_spacing_ref = 0.01;	/**< Initial reference particle spacing. */

Real rho0_f = 1000.0;						 /**< Reference density of fluid. */
Real rho0_v = 1.226;						 /**< Reference density of air. */
Real rho0_s = 7890;
Real gravity_g = 9.81;					 /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * 0.5); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
Real mu_water = 653.9e-6;
Real mu_air = 20.88e-6;

Real poisson = 0.27; 		/**< Poisson ratio.*/
Real Ae = 135.0e9; 			/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae;
std::string water = "./input/water.dat";
std::string tank_inner = "./input/tank_inner.dat";
std::string tank_outer = "./input/tank_outer.dat";


Real f = 1.0;
Real a = 0.08;

class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygonFromFile(water, ShapeBooleanOps::add);
		//multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for air block.
//----------------------------------------------------------------------
class VaporBlock : public MultiPolygonShape
{
public:
	explicit VaporBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygonFromFile(tank_inner, ShapeBooleanOps::add);
		multi_polygon_.addAPolygonFromFile(water, ShapeBooleanOps::sub);
		//multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Wall boundary shape definition.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{	
public:
	explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<MultiPolygonShape>(MultiPolygon(tank_outer), "OuterWall");
		subtract<MultiPolygonShape>(MultiPolygon(tank_inner), "InnerWall");
	}
};

class VariableGravity : public Gravity
{
	Real time_ = 0;
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g)) {};
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		time_ = GlobalStaticVariables::physical_time_;
		/*global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);*/
		global_acceleration_[0] = 0;
		return global_acceleration_;
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
	BoundingBox system_domain_bounds(Vec2d(-0.2, -0.090), Vec2d(0.2, 1.0));
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.setRunParticleRelaxation(false);
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.setReloadParticles(true);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f,mu_water);
	water_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<Vecd>("Acceleration");
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<Real>("Density");
	water_block.addBodyStateForRecording<Real>("DensitySummation");

	FluidBody vapor_block(sph_system, makeShared<VaporBlock>("VaporBody"));
	vapor_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_v, c_f,mu_air);
	vapor_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? vapor_block.generateParticles<ParticleGeneratorReload>(io_environment, vapor_block.getName())
		: vapor_block.generateParticles<ParticleGeneratorLattice>();
	

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	wall_boundary.defineComponentLevelSetShape("OuterWall");
	wall_boundary.defineComponentLevelSetShape("InnerWall");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
		: wall_boundary.generateParticles<ParticleGeneratorLattice>();
	//wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

	InnerRelation wall_boundary_inner(wall_boundary);
	ComplexRelation wall_boundary_complex(wall_boundary, RealBodyVector{&water_block, &vapor_block});
	InnerRelation water_inner(water_block);
	InnerRelation vapor_inner(vapor_block);


	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
		SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
		SimpleDynamics<RandomizeParticlePosition> random_vapor_particles(vapor_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_wall_boundary_to_vtp(io_environment, { &wall_boundary });
		BodyStatesRecordingToVtp write_water_to_vtp(io_environment, { &water_block });
		BodyStatesRecordingToVtp write_vapor_to_vtp(io_environment, { &vapor_block });
		/** Write the particle reload files. */
		ReloadParticleIO write_wall_boundary_particle_reload_files(io_environment, wall_boundary, "WallBoundary");
		ReloadParticleIO write_water_particle_reload_files(io_environment, water_block, "WaterBody");
		ReloadParticleIO write_vapor_particle_reload_files(io_environment, vapor_block, "VaporBody");
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
		relax_dynamics::RelaxationStepInner vapor_relaxation_step_inner(vapor_inner);
		relax_dynamics::RelaxationStepComplex wall_boundary_relaxation_step_complex(wall_boundary_complex, "OuterWall", true);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_wall_boundary_particles.exec(0.25);
		random_water_particles.exec(0.25);
		random_vapor_particles.exec(0.25);
		wall_boundary_relaxation_step_complex.SurfaceBounding().exec();
		water_relaxation_step_inner.SurfaceBounding().exec();
		vapor_relaxation_step_inner.SurfaceBounding().exec();
		write_wall_boundary_to_vtp.writeToFile(0);
		write_water_to_vtp.writeToFile(0);
		write_vapor_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			wall_boundary_relaxation_step_complex.exec();
			water_relaxation_step_inner.exec();
			vapor_relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
				write_water_to_vtp.writeToFile(ite_p);
				write_vapor_to_vtp.writeToFile(ite_p);
				write_wall_boundary_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_wall_boundary_particle_reload_files.writeToFile(0);
		write_water_particle_reload_files.writeToFile(0);
		write_vapor_particle_reload_files.writeToFile(0);

		return 0;
	}


	ComplexRelation water_vapor_complex(water_block, { &vapor_block });
	ContactRelation water_wall_contact(water_block, { &wall_boundary });
	ComplexRelation vapor_water_complex(vapor_block, { &water_block });
	ContactRelation vapor_wall_contact(vapor_block, { &wall_boundary });
	ContactRelation wall_boundary_water_contact(wall_boundary,RealBodyVector{&water_block,&vapor_block});
	//BodyRelationContact fluid_observer_contact(fluid_observer, RealBodyVector{ &water_block, &air_block });
	//BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	SimpleDynamics<NormalDirectionFromBodyShape> inner_normal_direction(wall_boundary);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> wall_boundary_corrected_configuration(wall_boundary_complex.getInnerRelation());
	//SimpleDynamics<ThermoVaporBodyInitialCondition> thermo_vapor_initial_condition(vapor_block);
	//SimpleDynamics<ThermoWaterBodyInitialCondition> thermo_water_initial_condition(water_block);
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_vapor_step(vapor_block, makeShared<VariableGravity>());
	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(water_block);
	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_vapor(vapor_block);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex>
		update_water_density_by_summation(water_wall_contact, water_vapor_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
		update_vapor_density_by_summation(vapor_wall_contact,vapor_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
		vapor_transport_correction(vapor_wall_contact,vapor_water_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_water(water_wall_contact, water_vapor_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_air(vapor_wall_contact, vapor_water_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_vapor_advection_time_step_size(vapor_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_vapor_time_step_size(vapor_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
		water_pressure_relaxation(water_wall_contact, water_vapor_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		water_density_relaxation(water_wall_contact, water_vapor_complex);
	/** Extend Pressure relaxation is used for air. */
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
		vapor_pressure_relaxation(vapor_wall_contact, vapor_water_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		vapor_density_relaxation(vapor_wall_contact, vapor_water_complex);
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RestartIO restart_io(io_environment, sph_system.real_bodies_);

	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	inner_normal_direction.exec();
	wall_boundary_corrected_configuration.exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
		water_block.updateCellLinkedList();
		vapor_block.updateCellLinkedList();
		water_vapor_complex.updateConfiguration();
		vapor_water_complex.updateConfiguration();
		water_wall_contact.updateConfiguration();
		vapor_wall_contact.updateConfiguration();

		//fluid_observer_contact.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 2.0;
	Real output_interval = 0.05;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	TickCount::interval_t interval_computing_time_step;
	TickCount::interval_t interval_computing_fluid_pressure_relaxation;
	TickCount::interval_t interval_updating_configuration;
	TickCount time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------

	//thermo_water_initial_condition.parallel_exec();
	//thermo_vapor_initial_condition.parallel_exec();
	body_states_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------

	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = TickCount::now();
			initialize_a_water_step.exec();
			initialize_a_vapor_step.exec();

			Real Dt_f = get_water_advection_time_step_size.exec();
			Real Dt_a = get_vapor_advection_time_step_size.exec();
			Real Dt = SMIN(Dt_f, Dt_a);

			update_water_density_by_summation.exec();
			update_vapor_density_by_summation.exec();
			//water_viscous_acceleration.parallel_exec();
			//vapor_viscous_acceleration.parallel_exec();
			vapor_transport_correction.exec();

			viscous_acceleration_air.exec();
			viscous_acceleration_water.exec();
			/*viscous_force_on_solid.parallel_exec();*/
		/*	wall_update_normal.parallel_exec();*/
			interval_computing_time_step += TickCount::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.exec();
				Real dt_a = get_vapor_time_step_size.exec();
				//Real dt_thermal_water = get_thermal_time_step_water.parallel_exec();
				//Real dt_thermal_air = get_thermal_time_step_vapor.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_a),Dt);

				water_pressure_relaxation.exec(dt);
				vapor_pressure_relaxation.exec(dt);
				/*fluid_force_on_solid_update.parallel_exec();*/
				water_density_relaxation.exec(dt);
				vapor_density_relaxation.exec(dt);

				size_t inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			vapor_block.updateCellLinkedListWithParticleSort(100);
			wall_boundary.updateCellLinkedList();
			vapor_water_complex.updateConfiguration();
			vapor_wall_contact.updateConfiguration();
			water_vapor_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();
			wall_boundary_water_contact.updateConfiguration();
			//fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += TickCount::now() - time_instance;
		}

		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}

	TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
		<< interval_computing_fluid_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";



	return 0;
};
