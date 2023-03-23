/**
 * @file 	two_phase_dambreak.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "two_phase_dambreak_static_confinement.h"
#include "sphinxsys.h"
#include "fluid_boundary_static_confinement.h"
using namespace SPH;

int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f);
	air_block.generateParticles<ParticleGeneratorLattice>();

	/*SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");*/

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_air_complex(water_block, {&air_block});
	//ContactRelation water_wall_contact(water_block, {&wall_boundary});
	ComplexRelation air_water_complex(air_block, {&water_block});
	//ContactRelation air_wall_contact(air_block, {&wall_boundary});
	ContactRelation fluid_observer_contact(fluid_observer, RealBodyVector{&water_block, &air_block});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));

	//SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(wall_boundary, "InnerWall");

	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, gravity_ptr);
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, gravity_ptr);
	/** Evaluation of density by summation approach. */
	//InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex>
	//	update_water_density_by_summation(water_wall_contact, water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_water_density_by_summation(water_air_complex.getInnerRelation());

	//InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
	//	update_air_density_by_summation(air_wall_contact, air_water_complex);
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_air_density_by_summation(air_water_complex);

	//InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
	//	air_transport_correction(air_wall_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>air_transport_correction(air_water_complex);

	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);

	/** Pressure relaxation for water by using position verlet time stepping. */
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
	//	water_pressure_relaxation(water_wall_contact, water_air_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemann> water_pressure_relaxation(water_air_complex);

	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
	//	water_density_relaxation(water_wall_contact, water_air_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemann> water_density_relaxation(water_air_complex);

	/** Extend Pressure relaxation is used for air. */
	//Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
	//	air_pressure_relaxation(air_wall_contact, air_water_complex, 2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemann>air_pressure_relaxation(air_water_complex);

	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
	//	air_density_relaxation(air_wall_contact, air_water_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemann>air_density_relaxation(air_water_complex);
	
	/** Confinement condition for wall and structure. */
	NearShapeSurface near_surface_water(water_block, makeShared<InnerWall>("InnerWall"));
	near_surface_water.level_set_shape_.writeLevelSet(io_environment);
	fluid_dynamics::StaticConfinementWithBounding confinement_condition_water(near_surface_water);

	NearShapeSurface near_surface_air(air_block, makeShared<InnerWall>("InnerWall"));
	fluid_dynamics::StaticConfinementWithPenalty confinement_condition_air(near_surface_air, c_f, 3.0);
	
	update_water_density_by_summation.post_processes_.push_back(&confinement_condition_water.density_summation_);
	water_pressure_relaxation.post_processes_.push_back(&confinement_condition_water.pressure_relaxation_);
	water_density_relaxation.post_processes_.push_back(&confinement_condition_water.density_relaxation_);
	//water_density_relaxation.post_processes_.push_back(&confinement_condition_water.surface_bounding_);

	update_air_density_by_summation.post_processes_.push_back(&confinement_condition_air.density_summation_);
	air_pressure_relaxation.post_processes_.push_back(&confinement_condition_air.extend_intergration_1st_half_);
	air_density_relaxation.post_processes_.push_back(&confinement_condition_air.density_relaxation_);
	//air_density_relaxation.post_processes_.push_back(&confinement_condition_air.surface_bounding_);
	air_transport_correction.post_processes_.push_back(&confinement_condition_air.transport_velocity_);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations 
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
	/** output the observed data from fluid body. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_pressure("Pressure", io_environment, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	//inner_normal_direction.parallel_exec();
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
 	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);
	write_recorded_pressure.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real end_time = 20.0;
	Real output_interval = 0.1;
	Real dt = 0.0;		  /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
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
			time_instance = tick_count::now();
			initialize_a_water_step.parallel_exec();
			initialize_a_air_step.parallel_exec();

			Real Dt_f = get_water_advection_time_step_size.parallel_exec();
			Real Dt_a = get_air_advection_time_step_size.parallel_exec();
			Real Dt = SMIN(Dt_f, Dt_a);

			update_water_density_by_summation.parallel_exec();
			update_air_density_by_summation.parallel_exec();

			air_transport_correction.parallel_exec();

			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.parallel_exec();
				Real dt_a = get_air_time_step_size.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_a), Dt);

				water_pressure_relaxation.parallel_exec(dt);
				air_pressure_relaxation.exec(dt);

				water_density_relaxation.parallel_exec(dt);
				air_density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
				{
					write_water_mechanical_energy.writeToFile(number_of_iterations);
					write_recorded_pressure.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_air_complex.updateConfiguration();
			//water_wall_contact.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_water_complex.updateConfiguration();
			//air_wall_contact.updateConfiguration();

			fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}

	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			  << interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
			  << interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	write_water_mechanical_energy.newResultTest();
	write_recorded_pressure.newResultTest();

	return 0;
}
