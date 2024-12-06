/**
 * @file 	heat_transfer_validation.cpp
 * @brief 	Heat Transfer in Slabs
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)

using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_r = 1000.0;					 /**< Reference density of right material. */
Real rho0_l = 1.226;						 /**< Reference density of left material. */
Real c_p_r = 4.179;
Real c_p_l = 1.012;
Real k_r = 0.620;
Real k_l = 0.0254;
Real diffusion_coff_r = k_r / (c_p_r * rho0_r);
Real diffusion_coff_l = k_l / (c_p_l * rho0_l);
Real dp = 0.0125;	/**< Initial reference particle spacing. */
Real initial_temperature_left = 0.0;
Real initial_temperature_rigth = 1.0;

Real c_ = 1.0;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric elements used in shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createOverallShape()
{
	std::vector<Vecd> over_all_shape;
	over_all_shape.push_back(Vecd(0.0, 0.0));
	over_all_shape.push_back(Vecd(0.0, 40.0 * dp));
	over_all_shape.push_back(Vecd(80.0 * dp, 40.0 * dp));
	over_all_shape.push_back(Vecd(80.0 * dp, 0.0));
	over_all_shape.push_back(Vecd(0.0, 0.0));

	return over_all_shape;
}

MultiPolygon createRightBlockShape()
{
    std::vector<Vecd> right_block_shape;
	right_block_shape.push_back(Vecd(40*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(right_block_shape, ShapeBooleanOps::add);

    return multi_polygon;
}


MultiPolygon createLeftBlockShape()
{
    std::vector<Vecd> left_block_shape;
	left_block_shape.push_back(Vecd(0.0, 0.0));
	left_block_shape.push_back(Vecd(0.0, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 0.0));
	left_block_shape.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(left_block_shape, ShapeBooleanOps::add);

    return multi_polygon;
}
class OverallBlock :public MultiPolygonShape
{
public:
	explicit OverallBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createOverallShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for right block.
//----------------------------------------------------------------------
//class RightBlock : public MultiPolygonShape
//{
//public:
//	explicit RightBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
//	{
//		multi_polygon_.addAPolygon(createRightBlockShape(), ShapeBooleanOps::add);
//	}
//};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for left block.
//----------------------------------------------------------------------



//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
template <class DynamicsIdentifier>
class LocalQuantityDefinition
    : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    LocalQuantityDefinition(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier){};
    virtual ~LocalQuantityDefinition(){};
};

class LocalDiffusivityDefinition : public LocalQuantityDefinition<BodyPartByParticle>
{
  public:
	  explicit LocalDiffusivityDefinition(BodyPartByParticle& body_part, Real local_diff, Real initial_phi)
		  : LocalQuantityDefinition<BodyPartByParticle>(body_part),
		  thermal_conductivity(particles_->getVariableDataByName<Real>("ThermalConductivity")),
		  phi(particles_->registerStateVariable<Real>("Phi")), initial_phi(initial_phi),
		  local_diff(local_diff)
	  {
		  particles_->addVariableToWrite<Real>("Phi");
	  };

    void update(size_t index_i, Real dt)
    {
        thermal_conductivity[index_i] = local_diff;
		phi[index_i] = initial_phi;
    };

  protected:
    Real *thermal_conductivity;
    Real local_diff;
	Real *phi;
	Real initial_phi;
};

//class DiffusionInitialCondition : public LocalDynamics
//{
//  public:
//    explicit DiffusionInitialCondition(SPHBody &sph_body)
//        : LocalDynamics(sph_body),
//          phi_(particles_->registerStateVariable<Real>("Phi")){};
//
//    void update(size_t index_i, Real dt)
//    {
//        phi_[index_i] = initial_temperature_left;
//    };
//
//  protected:
//    Real *phi_;
//};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
using ThermalRelaxationInner = DiffusionRelaxationRK2<
         DiffusionRelaxation<Inner<KernelGradientInner>, LocalIsotropicDiffusion>>;

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;

	size_t number_of_observation_points = 80;
		
		for(int i = 0;i< number_of_observation_points;i++)
		{
			observation_points.push_back(Vecd(i*dp,20*dp));
		}
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
	BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(80*dp, 40*dp));
	SPHSystem sph_system(system_domain_bounds, dp);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody diffusion_body(sph_system, makeShared<OverallBlock>("DiffusionBody"));
    LocalIsotropicDiffusion *overall_diffusion = diffusion_body.defineMaterial<LocalIsotropicDiffusion>("Phi", "Phi", diffusion_coff_l);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ContactRelation temperature_observer_contact(temperature_observer, {& diffusion_body });
	InnerRelation overall_inner(diffusion_body);
	// Define body regions
    BodyRegionByParticle right_body_particles(diffusion_body, makeShared<MultiPolygonShape>(createRightBlockShape()));
    BodyRegionByParticle left_body_particles(diffusion_body, makeShared<MultiPolygonShape>(createLeftBlockShape()));
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	 // Define diffusion coefficient
    SimpleDynamics<LocalDiffusivityDefinition> right_diffusivity(right_body_particles, diffusion_coff_r, 1.0);
    SimpleDynamics<LocalDiffusivityDefinition> left_diffusivity(left_body_particles, diffusion_coff_l, 0.0);

	ThermalRelaxationInner temperature_relaxation(overall_inner,overall_diffusion);

	GetDiffusionTimeStepSize get_time_step_size_right(diffusion_body, *overall_diffusion);
	//GetDiffusionTimeStepSize get_time_step_size_left(diffusion_body, left_diffusivity);
	//SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);

	/** Initialize particle acceleration. */
	//SimpleDynamics<ThermoRightBodyInitialCondition> thermo_right_initial_condition(right_block);
	//SimpleDynamics<ThermoLeftBodyInitialCondition> thermo_left_initial_condition(left_block);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	//RestartIO restart_io(io_environment, sph_system.real_bodies_);
	BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(diffusion_body, "Phi");
	ObservedQuantityRecording<Real> write_temperature("Phi", temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

    //setup_diffusion_initial_condition.exec();
   
    // thermal conductivity initialization
    right_diffusivity.exec();
    left_diffusivity.exec();

	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	int ite=0.0;
	Real end_time = 1.0;
	Real Output_Time = 0.1 * end_time;
	Real Observe_time =  0.05*Output_Time;
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
	//thermo_right_initial_condition.parallel_exec();
	//thermo_left_initial_condition.parallel_exec();
	write_real_body_states.writeToFile(0);
	write_temperature.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------

	while (physical_time < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < Output_Time)
		{
			time_instance = TickCount::now();
			interval_computing_time_step += TickCount::now() - time_instance;
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Observe_time)
			{
				Real dt_thermal_right = get_time_step_size_right.exec();
				//Real dt_thermal_left = get_time_step_size_left.exec();
				//dt = SMIN(dt_thermal_right, dt_thermal_left);
				dt = dt_thermal_right;
				temperature_relaxation.exec(dt);

				if (ite % 100== 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< physical_time << "	dt: "
						<< dt << "\n";
				}
				ite++;
				relaxation_time += dt;
				integration_time += dt;
				physical_time += dt;
			}

			interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;
			/** Update cell linked list and configuration. */
			diffusion_body.updateCellLinkedList();
			

			overall_inner.updateConfiguration();
			temperature_observer_contact.updateConfiguration();
			write_temperature.writeToFile();
			time_instance = TickCount::now();
			interval_updating_configuration += TickCount::now() - time_instance;
		}

		TickCount t2 = TickCount::now();
		write_real_body_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();
	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	return 0;
};
