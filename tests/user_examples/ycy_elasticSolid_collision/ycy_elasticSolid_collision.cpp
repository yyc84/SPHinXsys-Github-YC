/**
 * @file 	3d_elasticSolid_shell_collision.cpp
 * @brief 	This is a benchmark test of the 3D elastic solid->shell contact/impact formulations.
 * @details  We consider the collision of an elastic ball bouncing in a spherical shell box.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real PL = 4.0;				  /**< plate length. */
Real PW = 4.0;				  /**< plate width. */
Real SH = 4.0;					/**< system high. */
Real resolution_ref = 0.05;  /**< reference resolution. */
Real PH = 4* resolution_ref; /**< plate thickness for BCs. */
BoundingBox system_domain_bounds(Vec3d(0.0, 0.0, 0.0), Vec3d(PL, PW, SH));
Vec3d ball_center_1(2.0, 2.0, 2.0);
Real ball_radius = 0.5;
// observer location
StdVec<Vec3d> observation_location_1 = { ball_center_1 };
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real gravity_g = 1.0;
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
Real physical_viscosity = 10000.0;
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class Plate : public ComplexShape
{
public:
	explicit Plate(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vec3d halfsize_plate(0.5 * PL, 0.5 * PH, 0.5 * PW);
		Transform3d translation_plate(halfsize_plate);
		add<TransformShape<GeometricShapeBox>>(Transform3d(translation_plate), halfsize_plate);
	
	}
};


//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.setRunParticleRelaxation(false);
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.setReloadParticles(true);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------

	SolidBody ball(sph_system, makeShared<GeometricShapeBall>(ball_center_1, ball_radius, "BallBody"));
	ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	ball.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName())
		: ball.generateParticles<ParticleGeneratorLattice>();

	ObserverBody ball_observer(sph_system, "BallObserver");
	ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_1);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		InnerRelation ball_inner(ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation for ball.
		//----------------------------------------------------------------------
		SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(io_environment, ball);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		write_relaxed_particles.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.parallel_exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relaxed_particles.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	SolidBody plate(sph_system, makeShared<Plate>("plate"));
	plate.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
	//plate.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	plate.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation ball_inner(ball);
	SurfaceContactRelation ball_contact(ball, {&plate});
	ContactRelation ball_observer_contact(ball_observer, { &ball });
	////----------------------------------------------------------------------
	////	Define the main numerical methods used in the simulation.
	////	Note that there may be data dependence on the constructors of these methods.
	////----------------------------------------------------------------------
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g, 0.0));
	
	SimpleDynamics<TimeStepInitialization> ball_initialize_timestep(ball, gravity_ptr);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> ball_corrected_configuration(ball_inner);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size(ball);
	//** stress relaxation for the balls. */
	Dynamics1Level<solid_dynamics::KirchhoffIntegration1stHalf> ball_stress_relaxation_first_half(ball_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half(ball_inner);
	//** Algorithms for solid-solid contact. */
	InteractionDynamics<solid_dynamics::ContactDensitySummation, BodyPartByParticle> ball_update_contact_density(ball_contact);
	InteractionDynamics<solid_dynamics::ContactForceFromWall, BodyPartByParticle> ball_compute_solid_contact_forces(ball_contact);
	InteractionDynamics<solid_dynamics::DynamicContactForceWithWall, BodyPartByParticle> ball_plate_contact_force(ball_contact);
	
	

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	//BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	
	BodyStatesRecordingToVtp write_ball_state(io_environment, {ball});

	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	ObservedQuantityRecording<Vecd> ball_displacement_recording("Position", io_environment, ball_observer_contact);
	
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	ball_corrected_configuration.parallel_exec();
	//** Initial states output. */
	//write_ball_state.writeToFile(0);
	body_states_recording.writeToFile(0);
	ball_displacement_recording.writeToFile(0);

	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 10.0;
	Real end_time = T0;
	Real output_interval = 0.01 * T0;
	Real Dt = 0.1 * output_interval;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				ball_update_contact_density.parallel_exec();
				ball_compute_solid_contact_forces.parallel_exec();
				ball_plate_contact_force.parallel_exec();
				ball_stress_relaxation_first_half.parallel_exec(dt);
				//ball_friction.parallel_exec(dt);
				ball_stress_relaxation_second_half.parallel_exec(dt);

				ball.updateCellLinkedList();
				ball_contact.updateConfiguration();

				int num = ball_inner.base_particles_.pos_.size();
				int count_num = 0;
				for (size_t k = 0; k < ball_contact.contact_bodies_.size(); ++k)
				{
					for (int i = 0; i != num; ++i)
					{
						Neighborhood& neighborhood = ball_contact.contact_configuration_[k][i];
						if (neighborhood.current_size_ != 0)
						{
							count_num += 1;
						}
					}
				}
				
				std::string output_folder_ = "./output";
				std::string filefullpath = output_folder_ + "/" + "number_of_contact_particles" + ".dat";
				std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
				out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num << "  " << std::endl;


				ite++;
				Real dt_free = ball_get_time_step_size.parallel_exec();
				dt = dt_free;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				ball_displacement_recording.writeToFile(ite);
			}
		}
		tick_count t2 = tick_count::now();
		write_ball_state.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
