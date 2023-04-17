/**
 * @file 	dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#include "fluid_boundary_static_confinement.h"
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
class InnerWall : public ComplexShape
{	
public:
	explicit InnerWall(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<MultiPolygonShape>(MultiPolygon(tank_inner), "InnerWall");
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
	



	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------



	ComplexRelation water_vapor_complex(water_block, { &vapor_block });
	
	ComplexRelation vapor_water_complex(vapor_block, { &water_block });

	//BodyRelationContact fluid_observer_contact(fluid_observer, RealBodyVector{ &water_block, &air_block });
	//BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));

	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_vapor_step(vapor_block, makeShared<VariableGravity>());
	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(water_block);
	//GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_vapor(vapor_block);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner>
		update_water_density_by_summation(water_vapor_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
		update_vapor_density_by_summation(vapor_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
		vapor_transport_correction(vapor_water_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_vapor_advection_time_step_size(vapor_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_vapor_time_step_size(vapor_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalf>
		water_pressure_relaxation( water_vapor_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalf>
		water_density_relaxation(water_vapor_complex);
	/** Extend Pressure relaxation is used for air. */
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalf>
		vapor_pressure_relaxation(vapor_water_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalf>
		vapor_density_relaxation(vapor_water_complex);
	NearShapeSurface near_surface_water(water_block, makeShared<InnerWall>("InnerWall"));
	near_surface_water.level_set_shape_.writeLevelSet(io_environment);
	fluid_dynamics::StaticConfinementWithBounding confinement_condition_water(near_surface_water);

	NearShapeSurface near_surface_air(vapor_block, makeShared<InnerWall>("InnerWall"));
	fluid_dynamics::StaticConfinementWithPenalty confinement_condition_air(near_surface_air, c_f, 3.0);

	update_water_density_by_summation.post_processes_.push_back(&confinement_condition_water.density_summation_);
	water_pressure_relaxation.post_processes_.push_back(&confinement_condition_water.pressure_relaxation_);
	water_density_relaxation.post_processes_.push_back(&confinement_condition_water.density_relaxation_);
	//water_density_relaxation.post_processes_.push_back(&confinement_condition_water.surface_bounding_);

	update_vapor_density_by_summation.post_processes_.push_back(&confinement_condition_air.density_summation_);
	vapor_pressure_relaxation.post_processes_.push_back(&confinement_condition_air.extend_intergration_1st_half_);
	vapor_density_relaxation.post_processes_.push_back(&confinement_condition_air.density_relaxation_);
	//air_density_relaxation.post_processes_.push_back(&confinement_condition_air.surface_bounding_);
	vapor_transport_correction.post_processes_.push_back(&confinement_condition_air.transport_velocity_);
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
	/** output the observed data from fluid body. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);

	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	body_states_recording.writeToFile(0);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);

	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real end_time = 20.0;
	Real output_interval = 0.1;
	Real dt = 0.0;		  /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	TimeInterval interval_computing_time_step;
	TimeInterval interval_computing_pressure_relaxation;
	TimeInterval interval_updating_configuration;
	TickCount time_instance;
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

			vapor_transport_correction.exec();

			interval_computing_time_step += TickCount::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.exec();
				Real dt_a = get_vapor_time_step_size.exec();
				dt = SMIN(SMIN(dt_f, dt_a), Dt);

				water_pressure_relaxation.exec(dt);
				vapor_pressure_relaxation.exec(dt);

				water_density_relaxation.exec(dt);
				vapor_density_relaxation.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += TickCount::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
				{
					write_water_mechanical_energy.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_vapor_complex.updateConfiguration();
			//water_wall_contact.updateConfiguration();

			vapor_block.updateCellLinkedListWithParticleSort(100);
			vapor_water_complex.updateConfiguration();
			//air_wall_contact.updateConfiguration();

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
		<< interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";



	return 0;
};
