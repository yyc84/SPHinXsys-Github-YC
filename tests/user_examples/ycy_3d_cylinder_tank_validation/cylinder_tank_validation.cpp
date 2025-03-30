/**
 * @file 	sloshing_with_baffle.cpp
 * @brief 	This is a cylinder fall into static water.
 */
#include "sphinxsys.h"
#include "cylinder_tank_validation.h"
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
	// Tag for run particle relaxation for the initial body fitted distribution.
	system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	system.setReloadParticles(true);
	/* Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);


	/*
	@Brief creating body, materials and particles for the cylinder.
	*/
	SolidBody tank(system, makeShared<Tank>("Tank"));
	//tank.defineBodyLevelSetShape()->writeLevelSet(system);
	tank.defineMaterial<Solid>();
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? tank.generateParticles<BaseParticles, Reload>(tank.getName())
		: tank.generateParticles<BaseParticles, Lattice>();

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_water);
	//water_block.defineBodyLevelSetShape()->writeLevelSet(system);
	water_block.generateParticles<BaseParticles, Lattice>();

	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_a, c_f, mu_air);
	//air_block.defineBodyLevelSetShape()->writeLevelSet(system);
	air_block.generateParticles<BaseParticles, Lattice>();

	/*ObserverBody liquid_temperature_observer(system, "LiquidTemperatureObserver");
	liquid_temperature_observer.generateParticles<LiquidTemperatureObserverParticleGenerator>();

	ObserverBody gas_temperature_observer(system, "GasTemperatureObserver");
	gas_temperature_observer.generateParticles<GasTemperatureObserverParticleGenerator>();
	*/
	InnerRelation tank_inner(tank);
	InnerRelation water_inner(water_block);
	InnerRelation air_inner(air_block);

	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		using namespace relax_dynamics;
		SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
		//SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
		//SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_tank_to_vtp(tank );
		//BodyStatesRecordingToVtp write_water_to_vtp(water_block);
		//BodyStatesRecordingToVtp write_air_to_vtp(air_block );
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_particle_reload_files(tank);
		//ReloadParticleIO write_water_particle_reload_files(water_block);
		//ReloadParticleIO write_air_particle_reload_files(air_block);

		/** A  Physics relaxation step. */
		//relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
		//relax_dynamics::RelaxationStepInner air_relaxation_step_inner(air_inner);
		relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
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
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
				//write_water_to_vtp.writeToFile(ite_p);
				//write_air_to_vtp.writeToFile(ite_p);
				write_tank_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_tank_particle_reload_files.writeToFile(0);
		//write_water_particle_reload_files.writeToFile(0);
		//write_air_particle_reload_files.writeToFile(0);
		return 0;
	}


	ContactRelation water_tank_contact(water_block, { &tank });
	ContactRelation water_air_contact(water_block, {&air_block});
	ContactRelation air_water_contact(air_block, {&water_block});
	ContactRelation air_tank_contact(air_block, { &tank });
	
	//ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, { &water_block });
	//ContactRelation gas_temperature_observer_contact(gas_temperature_observer, { &air_block });
	
	//----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_tank_contact});
	
	
	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(system);
	
	RestartIO							restart_io(system);
	//ObservedQuantityRecording<Real> write_temperature_liquid("Phi", in_output, liquid_temperature_observer_contact);
	//ObservedQuantityRecording<Real> write_temperature_gas("Phi", in_output, gas_temperature_observer_contact);

	//InteractionDynamics<solid_dynamics::CorrectConfiguration> 		tank_corrected_configuration(tank_complex.inner_relation_);
	//Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    VariableGravity gravity(Vecd(0.0, 0.0, -gravity_g));
	
    SimpleDynamics<GravityForce> constant_gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce> constant_gravity_to_air(air_block, gravity);

	SimpleDynamics<NormalDirectionFromSubShapeAndOp> inner_normal_direction(tank, "InnerWall");
    InteractionDynamics<fluid_dynamics::BoundingFromWall> air_near_wall_bounding(air_tank_contact);
	

	/* Fluid dynamics */
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        water_pressure_relaxation(water_inner, water_air_contact, water_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        water_density_relaxation(water_inner, water_air_contact, water_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        air_pressure_relaxation(air_inner, air_water_contact, air_tank_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        air_density_relaxation(air_inner, air_water_contact, air_tank_contact);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface>
        update_water_density_by_summation(water_inner, water_tank_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_inner, air_water_contact, air_tank_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_inner, air_water_contact, air_tank_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_g);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);


	InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_acceleration_water(water_inner,  water_air_contact, water_tank_contact);
	InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_acceleration_air(air_inner, air_water_contact, air_tank_contact);

	DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec3d, FixedDampingRate>>>
        water_damping(0.2, ConstructorArgs(water_inner, "Velocity", mu_water), ConstructorArgs(water_air_contact, "Velocity", mu_water));
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec3d, FixedDampingRate>>>
        air_damping(0.2, ConstructorArgs(air_inner, "Velocity", mu_air), ConstructorArgs(air_water_contact, "Velocity", mu_air));
	//SimpleDynamics<ThermoAirBodyInitialCondition> thermo_air_initial_condition(air_block);
	//SimpleDynamics<ThermoWaterBodyInitialCondition> thermo_water_initial_condition(water_block);

	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(water_block);
	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_air(air_block);

	//ThermalRelaxationComplex thermal_relaxation_complex_water(water_air_complex);
	//ThermalRelaxationComplex thermal_relaxation_complex_air(air_water_complex);

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

	/*ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		water_average_temperature(in_output, water_block, "Phi");
	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		air_average_temperature(in_output, air_block, "Phi");
	ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>>
		air_rate_of_heat_transfer(in_output, air_block, "HeatFlux");
	ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>>
		water_rate_of_heat_transfer(in_output, water_block, "HeatFlux");*/

	//ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_air_total_mass(in_output, air_block, "MassiveMeasure");
	//ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_water_total_mass(in_output, water_block, "MassiveMeasure");
	
    ReducedQuantityRecording<TotalMechanicalEnergy> write_water_mechanical_energy(water_block, gravity);
	/**
	* @brief Pre-simulation.
	*/
	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	inner_normal_direction.exec();
	/** computing linear reproducing configuration for the tank. */
	//tank_corrected_configuration.parallel_exec();
	write_real_body_states.writeToFile(0);
	constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
	//thermo_water_initial_condition.parallel_exec();
	//thermo_air_initial_condition.parallel_exec();
	wave_probe_S1.writeToFile(0);
	wave_probe_S2.writeToFile(0);
	wave_probe_S3.writeToFile(0);
	write_real_body_states.addToWrite<Vecd>(tank, "NormalDirection"); 
	//water_average_temperature.writeToFile(0);
	//air_average_temperature.writeToFile(0);
	//air_rate_of_heat_transfer.writeToFile(0);
	//water_rate_of_heat_transfer.writeToFile(0);
	
	//compute_water_total_mass.writeToFile(0);
	//compute_air_total_mass.writeToFile(0);

	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		water_tank_contact.updateConfiguration();
		air_water_complex.updateConfiguration();
		air_tank_contact.updateConfiguration();
	}

	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 10.0;			/**< End time. */
	Real D_Time = 0.1;	/**< time stamps for output. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */


	/** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;


	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			time_instance = TickCount::now();
			/** outer loop for dual-time criteria time-stepping. */
			//initialize_a_water_step.parallel_exec();
			//initialize_a_air_step.parallel_exec();
			constant_gravity_to_water.exec();
			constant_gravity_to_air.exec();

			Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();

			Real Dt = SMIN(Dt_f, Dt_a);
			update_water_density_by_summation.exec();
            update_air_density_by_summation.exec();
            air_transport_correction.exec();
			viscous_acceleration_air.exec();
			viscous_acceleration_water.exec();
			air_near_wall_bounding.exec();

			interval_computing_time_step += TickCount::now() - time_instance;
			time_instance = TickCount::now();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
				//Real dt_thermal_water = get_thermal_time_step_water.parallel_exec();
				//Real dt_thermal_air = get_thermal_time_step_air.parallel_exec();
				//dt = SMIN(SMIN(dt_f, dt_thermal_water), SMIN(dt_thermal_air, dt_a), Dt);
				dt = SMIN(SMIN(dt_f, dt_a), Dt);
				/* Fluid pressure relaxation */
				
				/*if (GlobalStaticVariables::physical_time_ <= 5.0)
				{
                    water_damping.exec();
                    air_damping.exec();
				}*/
                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

				/* Fluid density relaxation */
				water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

				/*Thermal relaxation*/
				//thermal_relaxation_complex_air.parallel_exec(dt);
				//thermal_relaxation_complex_water.parallel_exec(dt);
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
			
			//liquid_temperature_observer_contact.updateConfiguration();
			//gas_temperature_observer_contact.updateConfiguration();
			//int num = water_air_complex.getInnerRelation().base_particles_.pos_.size();
			/*int count_num = 0;
			for (size_t k = 0; k < water_air_complex.contact_configuration_.size(); ++k)
			{
				for (int i = 0; i != num; ++i)
				{
					Neighborhood& neighborhood = water_air_complex.contact_configuration_[k][i];
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
				out_file_ << "\n";
				out_file_ << "Time " << "Num" << std::endl;
				out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num << "  " << std::endl;
			}*/

			number_of_iterations++;
			time_instance = TickCount::now();
			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedListWithParticleSort(100);
			//water_tank_contact.updateConfiguration();
			water_air_complex.updateConfiguration();
			air_block.updateCellLinkedListWithParticleSort(100);
			//air_tank_contact.updateConfiguration();
			air_water_complex.updateConfiguration();
			//write_real_body_states.writeToFile();
			interval_updating_configuration += TickCount::now() - time_instance;

			wave_probe_S1.writeToFile();
			wave_probe_S2.writeToFile();
            wave_probe_S3.writeToFile();
           
		}
		/** write run-time observation into file */
		TickCount t2 = TickCount::now();
		
		/*wave_probe_S1.writeToFile();
		wave_probe_S2.writeToFile();
		wave_probe_S3.writeToFile();*/
		//water_average_temperature.writeToFile();
		//air_average_temperature.writeToFile();
		//write_temperature_liquid.writeToFile();
		//write_temperature_gas.writeToFile();
		//air_rate_of_heat_transfer.writeToFile();
		//water_rate_of_heat_transfer.writeToFile();
		
		write_real_body_states.writeToFile();

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
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

	return 0;
}
