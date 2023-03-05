/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test to validate our anisotropic diffusion solver.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
using namespace SPH;   // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 2.0;
Real H = 0.4;
Real resolution_ref = H / 40.0;
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-4;
Real bias_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
std::array<std::string, 1> species_name_list{"Phi"};
//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
class DiffusionBlock : public MultiPolygonShape
{
public:
	explicit DiffusionBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> shape;
		shape.push_back(Vecd(0.0, 0.0));
		shape.push_back(Vecd(0.0, H));
		shape.push_back(Vecd(L, H));
		shape.push_back(Vecd(L, 0.0));
		shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
public:
	DiffusionMaterial() : DiffusionReaction<Solid>(species_name_list)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DiffusionInitialCondition
	: public DiffusionReactionInitialCondition<SolidParticles, Solid>
{
protected:
	size_t phi_;

public:
	explicit DiffusionInitialCondition(SPHBody &sph_body)
		: DiffusionReactionInitialCondition<SolidParticles, Solid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{

		if (pos_[index_i][0] >= 0.45 && pos_[index_i][0] <= 0.55)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		if (pos_[index_i][0] >= 1.0)
		{
			species_n_[phi_][index_i] = exp(-2500.0 * ((pos_[index_i][0] - 1.5) * (pos_[index_i][0] - 1.5)));
		}
	};
};
//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
class DiffusionBodyRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  RelaxationOfAllDiffusionSpeciesInner<SolidParticles, Solid>>
{
public:
	explicit DiffusionBodyRelaxation(InnerRelation &inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(inner_relation){};
	virtual ~DiffusionBodyRelaxation(){};
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class TemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit TemperatureObserverParticleGenerator(SPHBody &sph_body) 
	: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 11;
		Real range_of_measure = 0.9 * L;
		Real start_of_measure = 0.05 * L;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_of_observation_points - 1) + start_of_measure, 0.5 * H);
			positions_.push_back(point_coordinate);
		}
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody diffusion_body(sph_system, makeShared<DiffusionBlock>("DiffusionBlock"));
	diffusion_body.defineParticlesAndMaterial<DiffusionReactionParticles<SolidParticles, Solid>, DiffusionMaterial>();
	diffusion_body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Particle and body creation of fluid observers.
	//----------------------------------------------------------------------
	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
	temperature_observer.generateParticles<TemperatureObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation diffusion_body_inner_relation(diffusion_body);
	ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	DiffusionBodyRelaxation diffusion_relaxation(diffusion_body_inner_relation);
	SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> correct_configuration(diffusion_body_inner_relation);
	GetDiffusionTimeStepSize<SolidParticles, Solid> get_time_step_size(diffusion_body);
	PeriodicConditionUsingCellLinkedList periodic_condition_y(diffusion_body, diffusion_body.getBodyShapeBounds(), yAxis);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Real>>
		write_solid_temperature("Phi", io_environment, temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	correct_configuration.parallel_exec();
	setup_diffusion_initial_condition.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real end_time = T0;
	Real Output_Time = 0.1 * end_time;
	Real Observe_time = 0.1 * Output_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile();
	write_solid_temperature.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < Output_Time)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observe_time)
			{
				if (ite % 1 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}

				diffusion_relaxation.parallel_exec(dt);

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		write_solid_temperature.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_solid_temperature.newResultTest();

	return 0;
}
