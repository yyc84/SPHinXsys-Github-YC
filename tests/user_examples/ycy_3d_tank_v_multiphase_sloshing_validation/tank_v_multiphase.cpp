/**
 * @file 	waterentry.cpp
 * @brief 	This is a cylinder fall into static water.
 */
#include "sphinxsys.h"
#include "tank_case.h"
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
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for run particle relaxation for the initial body fitted distribution. */
    system.setRunParticleRelaxation(false);
	/** Tag for computation start with relaxed body fitted particles distribution. */
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
    tank.defineMaterial<Solid>();
	//tank.defineComponentLevelSetShape("OuterWall");
	//tank.defineComponentLevelSetShape("InnerWall");
	//tank.defineBodyLevelSetShape()->writeLevelSet(system);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
            ? tank.generateParticles<BaseParticles, Reload>(tank.getName())
            : tank.generateParticles<BaseParticles, Lattice>();

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_water);
	//water_block.defineBodyLevelSetShape()->writeLevelSet(system);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
            ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
            : water_block.generateParticles<BaseParticles, Lattice>();


	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
    air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_a, c_f, mu_air);
	//air_block.defineBodyLevelSetShape()->writeLevelSet(system);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
            ? air_block.generateParticles<BaseParticles, Reload>(air_block.getName())
            : air_block.generateParticles<BaseParticles, Lattice>();
	
	InnerRelation tank_inner(tank);
    InnerRelation water_inner(water_block);
    InnerRelation air_inner(air_block);

    ContactRelation tank_water_contact(tank, {&water_block});
    ContactRelation tank_air_contact(tank, {&air_block});
    ComplexRelation tank_complex(tank_inner, {&tank_water_contact, &tank_air_contact});
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		using namespace relax_dynamics;
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
		SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
		SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_tank_to_vtp(tank);
		BodyStatesRecordingToVtp write_water_to_vtp(water_block);
		BodyStatesRecordingToVtp write_air_to_vtp(air_block);
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_particle_reload_files(tank);
		ReloadParticleIO write_water_particle_reload_files(water_block);
		ReloadParticleIO write_air_particle_reload_files(air_block);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
		relax_dynamics::RelaxationStepInner air_relaxation_step_inner(air_inner);
        relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
		//relax_dynamics::RelaxationStepComplex tank_relaxation_step_complex(tank_complex, "OuterWall");
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_tank_particles.exec(0.25);
		random_water_particles.exec(0.25);
		random_air_particles.exec(0.25);
        //tank_relaxation_step_complex.SurfaceBounding().exec();
        water_relaxation_step_inner.SurfaceBounding().exec();
        air_relaxation_step_inner.SurfaceBounding().exec();
        tank_relaxation_step_inner.SurfaceBounding().exec();
		write_tank_to_vtp.writeToFile(0);
		write_water_to_vtp.writeToFile(0);
		write_air_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
            tank_relaxation_step_inner.exec();
			//tank_relaxation_step_complex.exec();
			water_relaxation_step_inner.exec();
			air_relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
				write_water_to_vtp.writeToFile(ite_p);
				write_air_to_vtp.writeToFile(ite_p);
				write_tank_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_tank_particle_reload_files.writeToFile(0);
		write_water_particle_reload_files.writeToFile(0);
		write_air_particle_reload_files.writeToFile(0);
		return 0;
	}


	/*ContactRelation water_block_contact(water_block,{ &tank });
	ContactRelation air_block_contact(air_block, { &tank });
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });*/

	ContactRelation water_tank_contact(water_block, {&tank});
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_tank_contact(air_block, {&tank});

	ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_tank_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_tank_contact});
	/*
	@Brief define simple data file input and outputs functions.
	*/
    BodyStatesRecordingToVtp write_real_body_states(system);
	
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

    //InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_acceleration_water(water_inner, water_air_contact, water_tank_contact);
    //InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> viscous_acceleration_air(air_inner, air_water_contact, air_tank_contact);



	//InteractionDynamics<solid_dynamics::CorrectConfiguration> 		tank_corrected_configuration(tank_complex.inner_relation_);
	//SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
	//SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	//SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, makeShared<VariableGravity>());
	/* Fluid dynamics */
	//InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(water_block_contact, water_air_complex.inner_relation_);
	//InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_block_contact,air_water_complex);
	//InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_block_contact,air_water_complex);
	//ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
	//ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> air_advection_time_step(air_block, U_g);
	//ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	//ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(air_block);
	//Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> fluid_pressure_relaxation(water_block_contact, water_air_complex.inner_relation_);
	//Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> fluid_density_relaxation( water_block_contact,water_air_complex.inner_relation_);
	//Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall> air_pressure_relaxation(air_block_contact, air_water_complex, 2.0);
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall> air_density_relaxation(air_block_contact,air_water_complex);
	
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
	/**
	 * @brief Pre-simulation.
	*/
	/** initialize cell linked lists for all bodies. */
    /** initialize cell linked lists for all bodies. */
    system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    system.initializeSystemConfigurations();
    /** computing surface normal direction for the tank. */
    inner_normal_direction.exec();
    /** computing linear reproducing configuration for the tank. */
    // tank_corrected_configuration.parallel_exec();
    write_real_body_states.writeToFile(0);
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
        
    // wave_probe_S1.writeToFile(0);
    // wave_probe_S2.writeToFile(0);
    // wave_probe_S3.writeToFile(0);
    write_real_body_states.addToWrite<Vecd>(tank, "NormalDirection");
       
	size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 20;
    Real End_Time = 22.0; /**< End time. */
    Real D_Time = 0.1;    /**< time stamps for output. */
    Real dt = 0.0;        /**< Default acoustic time step sizes for fluid. */

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
          
            Real Dt = SMIN(Dt_f, Dt_a);
            update_water_density_by_summation.exec();
            update_air_density_by_summation.exec();
            air_transport_correction.exec();
            //viscous_acceleration_air.exec();
            //viscous_acceleration_water.exec();
            air_near_wall_bounding.exec();

            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();

            Real relaxation_time = 0.0;

			while (relaxation_time < Dt)
			{
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;


				Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);
                /* Fluid pressure relaxation */
                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                /* Fluid density relaxation */
                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

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
			}
	
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
                write_water_mechanical_energy.writeToFile();
            }
            interval_updating_configuration += TickCount::now() - time_instance;
		}
		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
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
