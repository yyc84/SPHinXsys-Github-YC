/**
 * @file 	heat_transfer_validation.cpp
 * @brief 	Heat Transfer in Slabs
 */
#include "sphinxsys.h" //SPHinXsys Library.
//#include "heat_exchange_two_phase.h"
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
Real initial_temperature_right = 1.0;

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

std::vector<Vecd> createRightBlockShape()
{
    std::vector<Vecd> right_block_shape;
	right_block_shape.push_back(Vecd(40*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 0.0));

    return right_block_shape;
}


std::vector<Vecd> createLeftBlockShape()
{
    std::vector<Vecd> left_block_shape;
	left_block_shape.push_back(Vecd(0.0, 0.0));
	left_block_shape.push_back(Vecd(0.0, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 0.0));
	left_block_shape.push_back(Vecd(0.0, 0.0));

    return left_block_shape;
}


//----------------------------------------------------------------------
//	cases-dependent geometric shape for right block.
//----------------------------------------------------------------------
class RightBlock : public MultiPolygonShape
{
public:
	explicit RightBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createRightBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for left block.
//----------------------------------------------------------------------
class LeftBlock : public MultiPolygonShape
{
public:
	explicit LeftBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createLeftBlockShape(), ShapeBooleanOps::add);
	}
};


//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------

class LeftDiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit LeftDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>("Phi")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_left;
    };

  protected:
    Real *phi_;
};

class RightDiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit RightDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariable<Real>("Phi")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_right;
    };

  protected:
    Real *phi_;
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
using ThermalRelaxInner = DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>;
using ThermalRelaxContact = DiffusionRelaxation<Contact<KernelGradientContact>, IsotropicDiffusion, IsotropicDiffusion>;
//using ThermalRelaxContact = DiffusionRelaxation<HeatExchange<KernelGradientContact>, IsotropicDiffusion, IsotropicDiffusion>;
//using ThermalRelaxationComplex = TwoPhaseHeatExchangeBodyRelaxationComplex<HeatTransferDiffusion, HeatTransferDiffusion, KernelGradientInner, KernelGradientContact, TwoPhaseHeatExchange>;

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
    BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(80 * dp, 40 * dp));
    SPHSystem sph_system(system_domain_bounds, dp);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody right_body(sph_system, makeShared<RightBlock>("RightBody"));
    right_body.generateParticles<BaseParticles, Lattice>();

    FluidBody left_body(sph_system, makeShared<LeftBlock>("LeftBody"));
    right_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    ContactRelation temperature_observer_contact(temperature_observer, {&left_body, &right_body});
    InnerRelation left_inner(left_body);
    InnerRelation right_inner(right_body);
    ContactRelation left_body_contact(left_body, {&right_body});
    ContactRelation right_body_contact(right_body, {&left_body});

    ComplexRelation left_complex(left_inner, {&left_body_contact});
    ComplexRelation right_complex(right_inner, {&right_body_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // Define diffusion coefficient
    IsotropicDiffusion left_heat_diffusion("Phi", "Phi", diffusion_coff_l);
    IsotropicDiffusion right_heat_diffusion("Phi", "Phi", diffusion_coff_r);

    Dynamics1Level<ThermalRelaxInner> left_thermal_relax_inner(left_inner, &left_heat_diffusion);
    Dynamics1Level<ThermalRelaxInner> right_thermal_relax_inner(right_inner, &right_heat_diffusion);

    // Dynamics1Level<ThermalRelaxContact> left_thermal_relax_contact(left_body_contact, &left_heat_diffusion, {&right_heat_diffusion});
    // Dynamics1Level<ThermalRelaxContact> right_thermal_relax_contact(right_body_contact, &right_heat_diffusion, {&left_heat_diffusion});

    // ThermalRelaxContact thermal_diffusion(right_body_contact, &right_heat_diffusion, {&left_heat_diffusion});
    /*DiffusionRelaxation<Contact<KernelGradientContact>, IsotropicDiffusion> left_thermal_relax_contact(
            left_body_contact, &left_heat_diffusion);*/
    DiffusionRelaxation<Contact<KernelGradientContact>, IsotropicDiffusion, IsotropicDiffusion> left_thermal_relax_contact(
        left_body_contact, &left_heat_diffusion, &right_heat_diffusion);
    /*DiffusionRelaxation<TwoPhaseHeatExchange<KernelGradientContact>, IsotropicDiffusion, IsotropicDiffusion> left_thermal_relax_contact(
        left_body_contact, &left_heat_diffusion, {&right_heat_diffusion});*/

	/*ThermalRelaxationComplex thermal_relax_left_complex ( 
		ConstructorArgs(left_inner, &left_heat_diffusion),
        ConstructorArgs(left_body_contact, &left_heat_diffusion, &right_heat_diffusion));
	ThermalRelaxationComplex thermal_relax_right_complex ( 
		ConstructorArgs(right_inner, &right_heat_diffusion),
        ConstructorArgs(right_body_contact, &right_heat_diffusion, &left_heat_diffusion));*/

	SimpleDynamics<LeftDiffusionInitialCondition> left_diffusion_initial_condition(left_body);
	SimpleDynamics<RightDiffusionInitialCondition> right_diffusion_initial_condition(right_body);

	GetDiffusionTimeStepSize get_time_step_size_right(right_body, right_heat_diffusion);
	GetDiffusionTimeStepSize get_time_step_size_left(left_body, left_heat_diffusion);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	//RestartIO restart_io(io_environment, sph_system.real_bodies_);
	BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(right_body, "Phi");
	write_real_body_states.addToWrite<Real>(left_body, "Phi");
	ObservedQuantityRecording<Real> write_temperature("Phi", temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

   
    // thermal conductivity initialization
    left_diffusion_initial_condition.exec();
    right_diffusion_initial_condition.exec();

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
				Real dt_thermal_left = get_time_step_size_left.exec();
				dt = SMIN(dt_thermal_right, dt_thermal_left);
				left_thermal_relax_inner.exec(dt);
				right_thermal_relax_inner.exec(dt);
				//left_thermal_relax_contact.exec(dt);
				//right_thermal_relax_contact.exec(dt);
				//thermal_relax_left_complex.exec(dt);
				//thermal_relax_right_complex.exec(dt);

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
			left_body.updateCellLinkedList();
			right_body.updateCellLinkedList();

			left_complex.updateConfiguration();
			right_complex.updateConfiguration();

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
