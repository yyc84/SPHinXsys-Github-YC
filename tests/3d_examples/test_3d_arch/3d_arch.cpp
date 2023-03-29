/**
 * @file 	3d_arch.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider large deformation of a cylindrical thin structure.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 0.0975;								/** Radius of the inner boundary of the cylindrical thin structure. */
Real height = 0.02;									/** Height of the cylinder. */
Real thickness = 0.005;								/** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0; /** Radius of the mid surface. */
int particle_number = 10;							/** Particle number in the height direction. */
Real particle_spacing_ref = height / (Real)particle_number;
int particle_number_mid_surface = int(2.0 * radius_mid_surface * Pi * 215.0 / 360.0 / particle_spacing_ref);
int BWD = 1;								/** Width of the boundary layer measured by number of particles. */
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, 0.0, -radius - thickness),
								 Vec3d(radius + thickness, height, radius + thickness));
// Observer location
StdVec<Vecd> observation_location = {Vecd(0.0, height / 2.0, radius_mid_surface)};
/** For material properties of the solid. */
Real rho0_s = 7.800;			 /** Normalized density. */
Real Youngs_modulus = 210e6;	 /** Normalized Youngs Modulus. */
Real poisson = 0.3;				 /** Poisson ratio. */
Real physical_viscosity = 200.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.001;
Real gravitational_acceleration = -300000.0;

/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit CylinderParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the cylinder and boundary
		for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
		{
			for (int j = 0; j < particle_number; j++)
			{
				Real x = radius_mid_surface * cos(-17.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
				Real y = particle_spacing_ref * j + particle_spacing_ref * 0.5;
				Real z = radius_mid_surface * sin(-17.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
				Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
	}
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][2] < radius_mid_surface * sin(-17.5 / 180.0 * Pi))
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force ? current_time * global_acceleration_ / time_to_full_external_force : global_acceleration_;
	}
};

/**
 *  The main program
 */
int main(int ac, char *av[])
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	IOEnvironment io_environment(system);
	/** create a cylinder body with shell particles and linear elasticity. */
	SolidBody cylinder_body(system, makeShared<DefaultShape>("CylinderBody"));
	cylinder_body.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	cylinder_body.generateParticles<CylinderParticleGenerator>();

	/** Define Observer. */
	ObserverBody cylinder_observer(system, "CylinderObserver");
	cylinder_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation cylinder_body_inner(cylinder_body);
	ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder_body});

	/** Common particle dynamics. */
	SimpleDynamics<TimeStepInitialization> initialize_external_force(
		cylinder_body, makeShared<Gravity>(Vec3d(0.0, 0.0, gravitational_acceleration)));

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration(cylinder_body_inner);
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(cylinder_body);
	/** stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(cylinder_body_inner);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(cylinder_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(cylinder_body, "BoundaryGeometry");
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(boundary_geometry);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
		cylinder_position_damping(0.2, cylinder_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
		cylinder_rotation_damping(0.2, cylinder_body_inner, "AngularVelocity", physical_viscosity);
	/** Output */
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_cylinder_max_displacement("Position", io_environment, cylinder_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.exec();

	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_cylinder_max_displacement.writeToFile(0);

	/** Setup time stepping control parameters. */
	int ite = 0;
	Real end_time = 0.006;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			initialize_external_force.exec(dt);
			stress_relaxation_first_half.exec(dt);
			constrain_holder.exec(dt);
			cylinder_position_damping.exec(dt);
			cylinder_rotation_damping.exec(dt);
			constrain_holder.exec(dt);
			stress_relaxation_second_half.exec(dt);

			ite++;
			dt = computing_time_step_size.exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_cylinder_max_displacement.writeToFile(ite);
		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	if (system.generate_regression_data_)
	{
		write_cylinder_max_displacement.generateDataBase(0.005);
	}
	else
	{
		write_cylinder_max_displacement.newResultTest();
	}

	return 0;
}
