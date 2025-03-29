/**
 * @file 	sloshing_with_baffle.cpp
 * @brief 	This is a cylinder fall into static water.
 */
#include "sphinxsys.h"
#include "tank_case_heat.h"
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
        
	ObserverBody gas_temperature_observer(system, "GasTemperatureObserver");
    gas_temperature_observer.generateParticles<ObserverParticles>(GasTemperatureObserverParticle());

	ObserverBody liquid_temperature_observer(system, "LiquidTemperatureObserver");
    liquid_temperature_observer.generateParticles<ObserverParticles>(LiquidTemperatureObserverParticle());
	
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
	
	ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, { &water_block });
	ContactRelation gas_temperature_observer_contact(gas_temperature_observer, { &air_block });
	
	//----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_tank_contact});
	
	/*
	@Brief define simple data file input and outputs functions.
	*/

	//InteractionDynamics<solid_dynamics::CorrectConfiguration> 		tank_corrected_configuration(tank_complex.inner_relation_);
	VariableGravity gravity(Vecd(0.0, -gravity_g, 0.0));
	
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

	// Define diffusion coefficient
    HeatIsotropicDiffusion water_heat_diffusion("Phi", "Phi", k_water, rho0_f, c_p_water);
    HeatIsotropicDiffusion air_heat_diffusion("Phi", "Phi", k_air, rho0_a, c_p_air);

	Dynamics1Level<HeatExchangeComplex> water_heat_exchange_complex(water_inner, water_air_contact, &water_heat_diffusion, &air_heat_diffusion);
    Dynamics1Level<HeatExchangeComplex> air_heat_exchange_complex(air_inner, air_water_contact, &air_heat_diffusion, &water_heat_diffusion);

    SimpleDynamics<ThermoWaterBodyInitialCondition> water_diffusion_initial_condition(water_block);
    SimpleDynamics<ThermoAirBodyInitialCondition> air_diffusion_initial_condition(air_block);

	GetDiffusionTimeStepSize get_diffusion_time_step_size_water(water_block, water_heat_diffusion);
    GetDiffusionTimeStepSize get_diffusion_time_step_size_air(air_block, air_heat_diffusion);

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

	BodyStatesRecordingToVtp write_real_body_states(system);
    write_real_body_states.addToWrite<Real>(water_block, "Phi");
    write_real_body_states.addToWrite<Real>(air_block, "Phi");

    RestartIO restart_io(system);
    restart_io.addToWrite<Real>(water_block, "Phi");
    restart_io.addToWrite<Real>(air_block, "Phi");
    ObservedQuantityRecording<Real> write_temperature_liquid("Phi", liquid_temperature_observer_contact);
    ObservedQuantityRecording<Real> write_temperature_gas("Phi", gas_temperature_observer_contact);

	/*ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		water_average_temperature(in_output, water_block, "Phi");
	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		air_average_temperature(in_output, air_block, "Phi");*/
	
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_water_heat_flux_inner(water_block, "PhiFluxInner");
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_air_heat_flux_inner(air_block, "PhiFluxInner");
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_water_heat_flux_contact(water_block, "PhiFluxContact");
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_air_heat_flux_contact(air_block, "PhiFluxContact");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_water_heat_flux_wu(water_block, "PhiFluxWuContact");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_air_heat_flux_wu(air_block, "PhiFluxWuContact");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_water_heat_flux_inner_rate(water_block, "PhiFluxInnerChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_air_heat_flux_inner_rate(air_block, "PhiFluxInnerChangeRate");
    ReducedQuantityRecording<Average<QuantitySummation<Real, SPHBody>>> water_everage_temperature(water_block, "Phi");
    ReducedQuantityRecording<Average<QuantitySummation<Real, SPHBody>>> air_everage_temperature(air_block, "Phi");
    ReducedQuantityRecording<QuantityMax<Real, SPHBody>> water_max_temperature(water_block, "Phi");
    ReducedQuantityRecording<QuantityMax<Real, SPHBody>> air_max_temperature(air_block, "Phi");
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_water_heat_flux_contact_change_rate(water_block, "PhiChangeRate");
	/*ReducedQuantityRecording<QuantitySummation<Real,SPHBody>> compute_air_total_mass(air_block, "MassiveMeasure");
    ReducedQuantityRecording<QuantitySummation<Real,SPHBody>> compute_water_total_mass(water_block, "MassiveMeasure");*/
	
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
	water_diffusion_initial_condition.exec();
	air_diffusion_initial_condition.exec();
	//wave_probe_S1.writeToFile(0);
	//wave_probe_S2.writeToFile(0);
	//wave_probe_S3.writeToFile(0);
	write_real_body_states.addToWrite<Vecd>(tank, "NormalDirection"); 
	//water_average_temperature.writeToFile(0);
	//air_average_temperature.writeToFile(0);
	
	
	//compute_water_total_mass.writeToFile(0);
	//compute_air_total_mass.writeToFile(0);
    //write_water_heat_flux_inner.writeToFile(0);
   //write_air_heat_flux_inner.writeToFile(0);
    //write_water_heat_flux_contact.writeToFile(0);
    //write_air_heat_flux_contact.writeToFile(0);
    write_air_heat_flux_wu.writeToFile(0);
    write_water_heat_flux_wu.writeToFile(0);
    write_air_heat_flux_inner_rate.writeToFile(0);
    write_water_heat_flux_inner_rate.writeToFile(0);
    write_temperature_liquid.writeToFile(0);
    write_temperature_gas.writeToFile(0);
    //water_everage_temperature.writeToFile(0);
    //air_everage_temperature.writeToFile(0);
    if (system.RestartStep() != 0)
	{
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		//water_tank_contact.updateConfiguration();
		air_water_complex.updateConfiguration();
		//air_tank_contact.updateConfiguration();
	}

	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 20;
	Real End_Time = 22.0;			/**< End time. */
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
			constant_gravity_to_water.exec();
			constant_gravity_to_air.exec();

			Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real dt_water(0.0), dt_air(0.0);
            Real dt_f_thermal(0.0), dt_a_thermal(0.0);
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
				Real dt_thermal_water = get_diffusion_time_step_size_water.exec();
                Real dt_thermal_air = get_diffusion_time_step_size_air.exec();
				dt = SMIN(SMIN(dt_f, dt_thermal_water), SMIN(dt_thermal_air, dt_a), Dt);
				//dt = SMIN(SMIN(dt_f, dt_a), Dt);
                dt_water = dt_f;
                dt_air = dt_a;
                dt_f_thermal = dt_thermal_water;
                dt_a_thermal = dt_thermal_air;
				/* Fluid pressure relaxation */
                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

				/* Fluid density relaxation */
				water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

				/*Thermal relaxation*/
                if (GlobalStaticVariables::physical_time_>= 2.0)
                {
                    water_heat_exchange_complex.exec(dt);
                    air_heat_exchange_complex.exec(dt);
                }
                //water_heat_exchange_complex.exec(dt);
                //air_heat_exchange_complex.exec(dt);
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
                    << " Dt = " << Dt << " dt = " << dt << " Dt_water = " << Dt_f << " Dt_air = " << Dt_a
					<< " dt_f = " << dt_water << " dt_a = " << dt_air 
					<< " dt_f_thermal = " << dt_f_thermal << " dt_a_thermal = " << dt_a_thermal << "\n";

				/*if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);*/
				    
			}
			
			//liquid_temperature_observer_contact.updateConfiguration();
			//gas_temperature_observer_contact.updateConfiguration();
            int num_water = water_air_contact.getSPHBody().getBaseParticles().TotalRealParticles();
			int count_num_water = 0;
			for (size_t k = 0; k < water_air_contact.contact_configuration_.size(); ++k)
			{
                for (int i = 0; i != num_water; ++i)
				{
					Neighborhood& neighborhood = water_air_contact.contact_configuration_[k][i];
					if (neighborhood.current_size_ != 0)
					{
                        count_num_water += 1;
					}
				}
			}

			int num_air = air_water_contact.getSPHBody().getBaseParticles().TotalRealParticles();
            int count_num_air = 0;
            for (size_t k = 0; k < air_water_contact.contact_configuration_.size(); ++k)
            {
                for (int i = 0; i != num_air; ++i)
                {
                    Neighborhood &neighborhood = water_air_contact.contact_configuration_[k][i];
                    if (neighborhood.current_size_ != 0)
                    {
                        count_num_air += 1;
                    }
                }
            }
			
            std::string output_folder_ = "./output";
            std::string filefullpath = output_folder_ + "/" + "number_of_water-air_particles_has_contact" + ".dat";
            std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
            //out_file_ << "Time "<< "Num" << std::endl;        
			out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num_water << "  " << count_num_air << "  " << std::endl;

			number_of_iterations++;
			time_instance = TickCount::now();
			/** Update cell linked list and configuration. */
            water_block.updateCellLinkedList();
			air_block.updateCellLinkedList();
            water_air_complex.updateConfiguration();
            air_water_complex.updateConfiguration();

			if (GlobalStaticVariables::physical_time_ >= 2.0)
			{
                wave_probe_S1.writeToFile();
                wave_probe_S2.writeToFile();
                wave_probe_S3.writeToFile();
                write_water_heat_flux_inner.writeToFile();
                write_air_heat_flux_inner.writeToFile();
                write_water_heat_flux_contact.writeToFile();
                write_air_heat_flux_contact.writeToFile();
                write_air_heat_flux_wu.writeToFile();
                write_water_heat_flux_wu.writeToFile();
                write_air_heat_flux_inner_rate.writeToFile();
                write_water_heat_flux_inner_rate.writeToFile();
                write_temperature_liquid.writeToFile();
                write_temperature_gas.writeToFile();
                water_everage_temperature.writeToFile();
                air_everage_temperature.writeToFile();
                water_max_temperature.writeToFile();
                air_max_temperature.writeToFile();
                write_water_mechanical_energy.writeToFile();
                write_water_heat_flux_contact_change_rate.writeToFile();
			}
			
			//write_real_body_states.writeToFile();
			interval_updating_configuration += TickCount::now() - time_instance;
		}
		/** write run-time observation into file */
		TickCount t2 = TickCount::now();
		
		//wave_probe_S1.writeToFile();
		//wave_probe_S2.writeToFile();
		//wave_probe_S3.writeToFile();
        //compute_water_total_mass.writeToFile();
        //compute_air_total_mass.writeToFile();
        //write_water_heat_flux_inner.writeToFile();
        //write_air_heat_flux_inner.writeToFile();
        //write_water_heat_flux_contact.writeToFile();
        //write_air_heat_flux_contact.writeToFile();
        //write_air_heat_flux_wu.writeToFile();
        //write_water_heat_flux_wu.writeToFile();
        //write_air_heat_flux_inner_rate.writeToFile();
        //write_water_heat_flux_inner_rate.writeToFile();
		//water_average_temperature.writeToFile();
		//air_average_temperature.writeToFile();
		//write_temperature_liquid.writeToFile();
		//write_temperature_gas.writeToFile();
		
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
