/**
 * @file 	waterentry.cpp
 * @brief 	This is a cylinder fall into static water.
 */
#include "sphinxsys.h"
#include "tank_case.h"
#include "fluid_boundary_static_confinement.h"
 /*@brief Namespace cite here.
 */
using namespace SPH;



/*
Main program starts here.
*/
int main(int ac, char* av[])
{
	/* Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/* Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);


	/*
	@Brief creating body, materials and particles for the cylinder.
	*/

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<DiffusionWaterParticles , ThermoWaterBodyMaterial>();
	//water_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
	water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<Vecd>("Acceleration");
	water_block.addBodyStateForRecording<Vecd>("Acceleration");
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<Real>("Density");
	water_block.addBodyStateForRecording<Real>("DensitySummation");

	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<DiffusionAirParticles, ThermoAirBodyMaterial>();
	//air_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
	air_block.generateParticles<ParticleGeneratorLattice>();

	/*ObserverBody liquid_temperature_observer(system, "LiquidTemperatureObserver");
	liquid_temperature_observer.generateParticles<LiquidTemperatureObserverParticleGenerator>();

	ObserverBody gas_temperature_observer(system, "GasTemperatureObserver");
	gas_temperature_observer.generateParticles<GasTemperatureObserverParticleGenerator>();
	*/
	
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });


	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);
	//ObservedQuantityRecording<Real> write_temperature_liquid("Phi", in_output, liquid_temperature_observer_contact);
	//ObservedQuantityRecording<Real> write_temperature_gas("Phi", in_output, gas_temperature_observer_contact);

	//InteractionDynamics<solid_dynamics::CorrectConfiguration> 		tank_corrected_configuration(tank_complex.inner_relation_);
	SharedPtr<Gravity>gravity_ptr = makeShared<VariableGravity>();
	//SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, makeShared<VariableGravity>());
	/* Fluid dynamics */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> fluid_density_by_summation(water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> water_transport_correction(water_air_complex);

	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> air_advection_time_step(air_block, U_g);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(air_block);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> fluid_pressure_relaxation(water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> fluid_density_relaxation(water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemann> air_pressure_relaxation(air_water_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemann> air_density_relaxation(air_water_complex);
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_air_complex.getInnerRelation());
		/** Confinement condition for wall and structure. */
	NearShapeSurface near_surface_water(water_block, makeShared<WallAndStructure>("WallAndStructure"));
	near_surface_water.level_set_shape_.writeLevelSet(in_output);
	fluid_dynamics::StaticConfinementWithPenalty confinement_condition_water(near_surface_water, c_f, 2.0);

	NearShapeSurface near_surface_air(air_block, makeShared<WallAndStructure>("WallAndStructure"));
	fluid_dynamics::StaticConfinementWithPenalty confinement_condition_air(near_surface_air, c_f, 2.0);

	fluid_density_by_summation.post_processes_.push_back(&confinement_condition_water.density_summation_);
	fluid_pressure_relaxation.post_processes_.push_back(&confinement_condition_water.pressure_relaxation_);
	fluid_density_relaxation.post_processes_.push_back(&confinement_condition_water.density_relaxation_);
	/*fluid_density_relaxation.post_processes_.push_back(&confinement_condition_water.surface_bounding_);*/
	water_transport_correction.post_processes_.push_back(&confinement_condition_water.transport_velocity_);

	air_density_by_summation.post_processes_.push_back(&confinement_condition_air.density_summation_);
	air_pressure_relaxation.post_processes_.push_back(&confinement_condition_air.pressure_relaxation_);
	air_density_relaxation.post_processes_.push_back(&confinement_condition_air.density_relaxation_);
	air_density_relaxation.post_processes_.push_back(&confinement_condition_air.surface_bounding_);
	air_transport_correction.post_processes_.push_back(&confinement_condition_air.transport_velocity_);
	/*InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_water(water_block_contact, water_air_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_air(air_block_contact, air_water_complex);*/
	SimpleDynamics<ThermoAirBodyInitialCondition> thermo_air_initial_condition(air_block);
	SimpleDynamics<ThermoWaterBodyInitialCondition> thermo_water_initial_condition(water_block);

	TwoPhaseGetDiffusionTimeStepSize<DiffusionWaterParticles> get_thermal_time_step_water(water_block);
	TwoPhaseGetDiffusionTimeStepSize<DiffusionAirParticles> get_thermal_time_step_air(air_block);

	ThermalRelaxationComplex_1 thermal_relaxation_complex_water(water_air_complex);
	ThermalRelaxationComplex thermal_relaxation_complex_air(air_water_complex);

	BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("PorbeS1"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		probe_1(in_output, probe_s1);
	BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		probe_2(in_output, probe_s2);
	BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("PorbeS3"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		probe_3(in_output, probe_s3);

	//ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, DiffusionWaterParticles>>>
	//	water_average_temperature(in_output, water_block, "Phi");
	//ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, DiffusionAirParticles>>>
	//	air_average_temperature(in_output, air_block, "Phi");
	ReducedQuantityRecording<ReduceDynamics<QuantityMoment<Real>>>
		air_rate_of_heat_transfer(in_output, air_block, "HeatFlux");
	ReducedQuantityRecording<ReduceDynamics<QuantityMoment<Real>>>
		water_rate_of_heat_transfer(in_output, water_block, "HeatFlux");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_air_total_mass(in_output, air_block, "MassiveMeasure");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_water_total_mass(in_output, water_block, "MassiveMeasure");
	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */

	/** computing linear reproducing configuration for the tank. */
	//tank_corrected_configuration.parallel_exec();
	write_real_body_states.writeToFile(0);
	thermo_water_initial_condition.exec();
	thermo_air_initial_condition.exec();
	probe_1.writeToFile(0);
	probe_2.writeToFile(0);
	probe_3.writeToFile(0);
	/*water_average_temperature.writeToFile(0);
	air_average_temperature.writeToFile(0);*/
	air_rate_of_heat_transfer.writeToFile(0);
	water_rate_of_heat_transfer.writeToFile(0);
	compute_water_total_mass.writeToFile(0);
	compute_air_total_mass.writeToFile(0);
	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		air_water_complex.updateConfiguration();
		
	}

	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time =11.0;			/**< End time. */
	Real D_Time = 0.1;	/**< time stamps for output. */
	//Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */

	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			initialize_a_water_step.exec();
			initialize_a_air_step.exec();
			Real Dt_f = fluid_advection_time_step.exec();
			Real Dt_a = air_advection_time_step.exec();
			Real Dt = SMIN(Dt_f, Dt_a);
			fluid_density_by_summation.exec();
			air_density_by_summation.exec();
			air_transport_correction.exec();
			/*viscous_acceleration_air.exec();
			viscous_acceleration_water.exec();*/

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = fluid_acoustic_time_step.exec();
				Real dt_a = air_acoustic_time_step.exec();
				Real dt_thermal_water = get_thermal_time_step_water.exec();
				Real dt_thermal_air = get_thermal_time_step_air.exec();
				dt = SMIN(SMIN(dt_f, dt_thermal_water), SMIN(dt_thermal_air, dt_a), Dt);
				/* Fluid pressure relaxation */
				fluid_pressure_relaxation.exec(dt);
				air_pressure_relaxation.exec(dt);
				/* Fluid density relaxation */
				fluid_density_relaxation.exec(dt);
				air_density_relaxation.exec(dt);
				thermal_relaxation_complex_air.exec(dt);
				thermal_relaxation_complex_water.exec(dt);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}

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

			/** Update cell linked list and configuration. */

			water_block.updateCellLinkedListWithParticleSort(100);
			water_air_complex.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_water_complex.updateConfiguration();
			air_rate_of_heat_transfer.writeToFile();
			water_rate_of_heat_transfer.writeToFile();
			/*liquid_temperature_observer_contact.updateConfiguration();
			gas_temperature_observer_contact.updateConfiguration();*/
			int num = air_water_complex.getInnerRelation().base_particles_.pos_.size();
			int count_num = 0;
			for (size_t k = 0; k < air_water_complex.contact_configuration_.size(); ++k)
			{
				for (int i = 0; i != num; ++i)
				{
					Neighborhood& neighborhood = air_water_complex.contact_configuration_[k][i];
					if (neighborhood.current_size_ != 0)
					{
						count_num += 1;
					}
				}
			}
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::string output_folder_ = "./output";
				std::string filefullpath = output_folder_ + "/" + "number_of_contact_particles" + ".dat";
				std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
				out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num << "  " << std::endl;
			}
		}
		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
		compute_vorticity.exec();
		write_real_body_states.writeToFile();
		probe_1.writeToFile();
		probe_2.writeToFile();
		probe_3.writeToFile();
		/*water_average_temperature.writeToFile();
		air_average_temperature.writeToFile();*/
		//write_temperature_liquid.writeToFile();
		//write_temperature_gas.writeToFile();



		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
