/**
 * @file 	heat_transfer_validation.cpp
 * @brief 	Heat Transfer in Slabs
 */
#include "sphinxsys.h" //SPHinXsys Library.
//#include "heat_diffusion_dynamics.h"
//#include "heat_diffusion_dynamics.hpp"
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

class LeftDiffusionInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit LeftDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
		heat_flux_inner_(*particles_->registerSharedVariable<Real>("PhiFluxInner")),
          heat_flux_contact_(*particles_->registerSharedVariable<Real>("PhiFluxContact")),
          heat_flux_wu_contact_(*particles_->registerSharedVariable<Real>("PhiFluxWuContact")) {};
    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_left;
    };

  protected:
    StdLargeVec<Real> &phi_;
	StdLargeVec<Real> &heat_flux_inner_;
    StdLargeVec<Real> &heat_flux_contact_;
    StdLargeVec<Real> &heat_flux_wu_contact_;
    
};

class RightDiffusionInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit RightDiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
		heat_flux_inner_(*particles_->registerSharedVariable<Real>("PhiFluxInner")),
          heat_flux_contact_(*particles_->registerSharedVariable<Real>("PhiFluxContact")),
          heat_flux_wu_contact_(*particles_->registerSharedVariable<Real>("PhiFluxWuContact")) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_right;
    };

  protected:
    StdLargeVec<Real> &phi_;
	StdLargeVec<Real> &heat_flux_inner_;
    StdLargeVec<Real> &heat_flux_contact_;
    StdLargeVec<Real> &heat_flux_wu_contact_;
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
//using ThermalRelaxInner = HeatDiffusionRelaxation<HeatInner<KernelGradientInner>, HeatIsotropicDiffusion>;

using ThermalDiffusionContact = DiffusionRelaxation<Contact<KernelGradientContact>, HeatIsotropicDiffusion,HeatIsotropicDiffusion>;
using ThermalDiffusionInner = DiffusionRelaxation<HeatInner<KernelGradientInner>, HeatIsotropicDiffusion>;

using ThermalDiffusionComplex = HeatExchangeDiffusionComplex<KernelGradientInner, KernelGradientContact, HeatIsotropicDiffusion, HeatIsotropicDiffusion>;
//using ThermalRelaxationComplex = HeatDiffusionBodyRelaxationComplex<HeatIsotropicDiffusion,HeatIsotropicDiffusion, KernelGradientInner, KernelGradientContact, Contact>;

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;

	size_t number_of_observation_points = 160;
		
		for(int i = 0;i< number_of_observation_points;i++)
		{
			observation_points.push_back(Vecd(i*0.5*dp,20*dp));
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
    left_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    ContactRelation temperature_observer_contact(temperature_observer, {&left_body, &right_body});
    InnerRelation right_inner(right_body);
	InnerRelation left_inner(left_body);
    
    ContactRelation left_body_contact(left_body, {&right_body});
    ContactRelation right_body_contact(right_body, {&left_body});

    ComplexRelation left_complex(left_inner, {&left_body_contact});
    ComplexRelation right_complex(right_inner, {&right_body_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // Define diffusion coefficient
    HeatIsotropicDiffusion left_heat_diffusion("Phi", "Phi", k_l, rho0_l, c_p_l);
    HeatIsotropicDiffusion right_heat_diffusion("Phi", "Phi", k_r, rho0_r, c_p_r);
    IsotropicDiffusion right_diffusion("Phi", "Phi", k_r);

    
    //ThermalDiffusionInner right_thermal_diffusion_inner(right_inner, &right_heat_diffusion);
    //ThermalDiffusionContact right_thermal_diffusion_contact(right_body_contact, &right_heat_diffusion, &left_heat_diffusion);
    //ThermalDiffusionComplex right_thermal_diffusion_complex(right_inner, right_body_contact, &right_heat_diffusion, &left_heat_diffusion);

    //ThermalDiffusionInner left_thermal_diffusion_inner(left_inner, &left_heat_diffusion);
    //ThermalDiffusionContact left_thermal_diffusion_contact(left_body_contact, &left_heat_diffusion, &right_heat_diffusion);

	//Dynamics1Level<ThermalDiffusionInner> left_thermal_relax_inner(left_inner, &left_heat_diffusion);
	//Dynamics1Level<ThermalDiffusionContact> left_thermal_relax_contact(left_body_contact, &left_heat_diffusion, &right_heat_diffusion);
    
    Dynamics1Level<ThermalDiffusionComplex,SequencedPolicy> left_thermal_relax_complex(left_inner, left_body_contact, &left_heat_diffusion, &right_heat_diffusion);
    Dynamics1Level<ThermalDiffusionComplex,SequencedPolicy> right_thermal_relax_complex(right_inner, right_body_contact, &right_heat_diffusion, &left_heat_diffusion);

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
	ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_change_rate_total(right_body, "PhiChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_change_rate_total(left_body, "PhiChangeRate");
	ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_change_rate_inner(right_body, "PhiFluxInnerChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_change_rate_inner(left_body, "PhiFluxInnerChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_inner(right_body, "PhiFluxInner");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_inner(left_body, "PhiFluxInner");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_change_rate_contact(right_body, "PhiFluxContactChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_change_rate_contact(left_body, "PhiFluxContactChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_contact(right_body, "PhiFluxContact");
    ReducedQuantityRecording<QuantitySummation<Real, SPHBody>> write_left_heat_flux_contact(left_body, "PhiFluxContact");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_change_rate_contact_wu(right_body, "PhiFluxWuContactChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_change_rate_contact_wu(left_body, "PhiFluxWuContactChangeRate");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_right_heat_flux_contact_wu(right_body, "PhiFluxWuContact");
    ReducedQuantityRecording<QuantityMoment<Real, SPHBody>> write_left_heat_flux_contact_wu(left_body, "PhiFluxWuContact");
    
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
	//Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
	write_right_heat_flux_inner.writeToFile(0);
	write_left_heat_flux_inner.writeToFile(0);
    write_right_heat_flux_contact.writeToFile(0);
    write_left_heat_flux_contact.writeToFile(0);
    write_right_heat_flux_contact_wu.writeToFile(0);
    write_left_heat_flux_contact_wu.writeToFile(0);
    write_right_heat_flux_change_rate_total.writeToFile(0);
    write_left_heat_flux_change_rate_total.writeToFile(0);
    write_right_heat_flux_change_rate_inner.writeToFile(0);
    write_left_heat_flux_change_rate_inner.writeToFile(0);
    write_right_heat_flux_change_rate_contact.writeToFile(0);
    write_left_heat_flux_change_rate_contact.writeToFile(0);
    write_right_heat_flux_change_rate_contact_wu.writeToFile(0);
    write_left_heat_flux_change_rate_contact_wu.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------

	while (GlobalStaticVariables::physical_time_ < end_time)
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
				dt = 0.2*SMIN(dt_thermal_right, dt_thermal_left);
                right_thermal_relax_complex.exec(dt);
                left_thermal_relax_complex.exec(dt);

				if (ite % 100== 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				ite++;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				write_temperature.writeToFile();
			}

			interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;
			/** Update cell linked list and configuration. */
			left_body.updateCellLinkedList();
			right_body.updateCellLinkedList();

			left_complex.updateConfiguration();
			right_complex.updateConfiguration();

			temperature_observer_contact.updateConfiguration();
			//write_temperature.writeToFile();
			time_instance = TickCount::now();
			interval_updating_configuration += TickCount::now() - time_instance;
		}

		TickCount t2 = TickCount::now();
		write_real_body_states.writeToFile();
		write_right_heat_flux_inner.writeToFile();
		write_left_heat_flux_inner.writeToFile();
		write_right_heat_flux_contact.writeToFile();
		write_left_heat_flux_contact.writeToFile();
        write_right_heat_flux_contact_wu.writeToFile();
        write_left_heat_flux_contact_wu.writeToFile();
        write_right_heat_flux_change_rate_total.writeToFile();
        write_left_heat_flux_change_rate_total.writeToFile();
        write_right_heat_flux_change_rate_inner.writeToFile();
        write_left_heat_flux_change_rate_inner.writeToFile();
        write_right_heat_flux_change_rate_contact.writeToFile();
        write_left_heat_flux_change_rate_contact.writeToFile();
        write_right_heat_flux_change_rate_contact_wu.writeToFile();
        write_left_heat_flux_change_rate_contact_wu.writeToFile();
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
