/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;						  /**< Reference length. */
Real DH = 3.0;						  /**< Reference and the height of main channel. */
Real DL1 = 2.2 * DL;				  /**< The length of the main channel. */
Real DL_i = 0.2 * DL;                 /**< The width of injection */
Real DH_i = 0.5 * DH;                 /**< The high of injection */
Real resolution_ref = 0.15;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real rho0_a = 0.001;
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(1.0, DH / (2.0 * (DL - DL1)));
Real Re = 100.0;					/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
Real mu_a = rho0_a * U_f * DH / Re; 
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the water block for single phase flow. */
//std::vector<Vecd> water_block_shape{
//	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 0.0),
//	Vecd(DL + DL_i, 0.0), Vecd(DL + DL_i, -DH_i), Vecd(DL, -DH_i), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0) };
///** the water block in for two phase flow. */
//std::vector<Vecd> water_block_shape_1{
//	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0)};
///** the air block in T shape polygon. */
//std::vector<Vecd> air_block_shape{
//	Vecd(DL, -DH_i), Vecd(DL, 0.0), Vecd(DL + DL_i, 0.0), Vecd(DL + DL_i, -DH_i), };
///** the outer wall polygon. */
//std::vector<Vecd> outer_wall_shape{
//	Vecd(-DL_sponge - BW, -BW), Vecd(-DL_sponge - BW, DH + BW), Vecd(DL1 + BW, DH + BW), Vecd(DL1 + BW, -BW), Vecd(DL + DL_i + BW, -BW), 
//	Vecd(DL + DL_i + BW, -DH_i - BW), Vecd(DL - BW, -DH_i - BW), Vecd(DL - BW, -BW), Vecd(-DL_sponge - BW, -BW)};
///** the inner wall polygon. */
//std::vector<Vecd> inner_wall_shape{
//	Vecd(-DL_sponge - BW, 0.0), Vecd(-DL_sponge - BW, DH), Vecd(DL1+BW, DH), Vecd(DL1+BW, 0.0), Vecd(DL+DL_i, 0.0), 
//	Vecd(DL+DL_i, -DH_i - BW), Vecd(DL, -DH_i - BW), Vecd(DL, 0.0), Vecd(-DL_sponge - BW, 0.0)};

/*for test== injection is on upside*/
//std::vector<Vecd> water_block_shape{
//	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, DH + DH_i), Vecd(DL + DL_i, DH + DH_i), Vecd(DL + DL_i, DH),
//	Vecd(DL1, DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0) };
///** the water block in for two phase flow. */
//std::vector<Vecd> water_block_shape_1{
//	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0) };
///** the air block in T shape polygon. */
//std::vector<Vecd> air_block_shape{
//	Vecd(DL, DH), Vecd(DL, DH + DH_i), Vecd(DL + DL_i, DH + DH_i), Vecd(DL + DL_i, DH) };
///** the outer wall polygon. */
//std::vector<Vecd> outer_wall_shape{
//	Vecd(-DL_sponge - BW, -BW), Vecd(-DL_sponge - BW, DH + BW), Vecd(DL-BW, DH + BW), Vecd(DL-BW, DH + DH_i + BW), 
//	Vecd(DL + DL_i + BW, DH + DH_i + BW), Vecd(DL + DL_i + BW, DH + BW),
//	Vecd(DL1 + BW, DH + BW), Vecd(DL1 + BW, -BW), Vecd(-DL_sponge - BW, -BW) };
///** the inner wall polygon. */
//std::vector<Vecd> inner_wall_shape{
//	Vecd(-DL_sponge - BW, 0.0), Vecd(-DL_sponge - BW, DH), Vecd(DL, DH), Vecd(DL, DH + DH_i + BW), Vecd(DL + DL_i, DH + DH_i + BW), Vecd(DL + DL_i, DH),
//	Vecd(DL1 + BW, DH), Vecd(DL1 + BW, 0.0),  Vecd(-DL_sponge - BW, 0.0) };

/*for test two injection from left*/
std::vector<Vecd> water_block_shape{
	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, 1.5 * DH + 2 * DH_i), Vecd(DL1, 1.5 * DH + 2 * DH_i), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0) };
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
	Vecd(-DL_sponge - BW, -BW), Vecd(-DL_sponge - BW, 1.5 * DH + 2 * DH_i + BW), Vecd(DL1 + BW, 1.5 * DH + 2 * DH_i + BW), Vecd(DL1 + BW, -BW), Vecd(-DL_sponge - BW, -BW) };
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, 1.5 * DH + 2 * DH_i), Vecd(DL1 + BW, 1.5 * DH + 2 * DH_i), Vecd(DL1 + BW, 0.0),  Vecd(-DL_sponge, 0.0) };
std::vector<Vecd> injection1{
	Vecd(-DL_sponge - BW, 0.5 * DH), Vecd(-DL_sponge - BW, 0.5 * DH + DH_i), Vecd(-DL_sponge, 0.5 * DH + DH_i), Vecd(-DL_sponge, 0.5 * DH),Vecd(-DL_sponge - BW, 0.5 * DH)};
std::vector<Vecd> injection2{
	Vecd(-DL_sponge - BW, DH + DH_i), Vecd(-DL_sponge - BW, DH + 2 * DH_i), Vecd(-DL_sponge, DH + 2 * DH_i), Vecd(-DL_sponge, DH + DH_i), Vecd(-DL_sponge - BW, DH + DH_i) };
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(injection1, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(injection2, ShapeBooleanOps::add);
	}
};

class AirBlock : public MultiPolygonShape
{
public:
	explicit AirBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(injection1, ShapeBooleanOps::add);
	}
};

class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		multi_polygon_.addAPolygon(injection1, ShapeBooleanOps::sub);
		multi_polygon_.addAPolygon(injection2, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
	Real u_ref_, t_ref_;
	AlignedBoxShape &aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity(BoundaryConditionType &boundary_condition)
		: u_ref_(U_f), t_ref_(2.0),
		  aligned_box_(boundary_condition.getAlignedBox()),
		  halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd &position, Vecd &velocity)
	{
		Vecd target_velocity = velocity;
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		if (aligned_box_.checkInBounds(0, position))
		{
			target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
		}
		return target_velocity;
	}
};

/*for test*/
struct InflowVelocity1
{
	Real u_ref_, t_ref_;
	AlignedBoxShape& aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity1(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(2.0),
		aligned_box_(boundary_condition.getAlignedBox()),
		halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = velocity;
		//Real run_time = GlobalStaticVariables::physical_time_;
		//Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		if (aligned_box_.checkInBounds(0, position))
		{
			target_velocity[1] = 1.0;
		}
		return target_velocity;
	}
};
/*for pipe injection*/
struct InflowVelocity2
{
	Real u_ref_, t_ref_;
	AlignedBoxShape& aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity2(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(2.0),
		aligned_box_(boundary_condition.getAlignedBox()),
		halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = velocity;
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		if (aligned_box_.checkInBounds(0, position))
		{
			target_velocity[1] = 1.5 * u_ave * (1.0 - position[0] * position[0] / halfsize_[0] / halfsize_[0]);
		}
		return target_velocity;
	}
};

/*for test*/
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
public:
	InletInflowCondition(BodyAlignedBoxByParticle& aligned_box_part)
		: EmitterInflowCondition(aligned_box_part) {}

protected:
	virtual Vecd getTargetVelocity(Vecd& position, Vecd& velocity) override
	{
		return Vec2d(2.0, 0.0);
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -DH_i - BW), Vec2d(DL1 + BW, 2 * (DH + DH_i + BW)));
	SPHSystem system(system_domain_bounds, resolution_ref);
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
	water_block.generateParticles<ParticleGeneratorLattice>();

	/* for two phase flow */
	/*FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
	air_block.generateParticles<ParticleGeneratorLattice>();*/

	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

	//BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);
	///write_body_states.writeToFile();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	/*for single phase flow*/
	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex_relation(water_block_inner, {&wall_boundary });

	/*for two phase flow*/
	/*ComplexRelation water_air_complex(water_block, { &air_block });
	ContactRelation water_wall_contact(water_block, { &wall_boundary });
	ComplexRelation air_water_complex(air_block, { &water_block });
	ContactRelation air_wall_contact(air_block, { &wall_boundary });*/


	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/*for single phase flow*/
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation_water(water_block_complex_relation);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation_water(water_block_complex_relation);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_water(water_block_complex_relation);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction_water(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
		inlet_outlet_surface_particle_indicator_water(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation_water(water_block_complex_relation);
	water_block.addBodyStateForRecording<Real>("Pressure");		   // output for debug
	water_block.addBodyStateForRecording<int>("SurfaceIndicator"); // output for debug

	/*for single phase flow*/
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step_water(water_block);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size_water(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size_water(water_block);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

	/*for two phase flow*/
	//Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation_air(air_block_complex_relation);
	//Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation_air(air_block_complex_relation);
	//InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_air(air_block_complex_relation);
	//InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction_air(air_block_complex_relation);
	//InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
	//	inlet_outlet_surface_particle_indicator_air(water_block_complex_relation);
	//InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation_air(air_block_complex_relation);
	//air_block.addBodyStateForRecording<Real>("Pressure");		   // output for debug
	//air_block.addBodyStateForRecording<int>("SurfaceIndicator"); // output for debug

	//SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	////SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(wall_boundary, "InnerWall");
	//SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step_water(water_block);
	//SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step_air(air_block);
	///** Time step size without considering sound wave speed. */
	//ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size_water(water_block, U_f);
	//ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size_air(air_block, U_f);
	///** Time step size with considering sound wave speed. */
	//ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size_water(water_block);
	//ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size_air(air_block);
	///** Evaluation of density by summation approach. */
	//InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex>
	//	update_density_by_summation_water(water_wall_contact, water_air_complex.getInnerRelation());
	//InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex>
	//	update_density_by_summation_air(air_wall_contact, air_water_complex.getInnerRelation());
	//InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
	//	transport_velocity_correction_water(water_wall_contact, water_air_complex);
	//InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
	//	transport_velocity_correction_air(air_wall_contact, air_water_complex);
	///*surface identification*/
	//InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
	//	inlet_outlet_surface_particle_indicator_water(water_air_complex.getInnerRelation(), water_wall_contact);
	//InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
	//	inlet_outlet_surface_particle_indicator_air(air_water_complex.getInnerRelation(), air_wall_contact);
	///*Viscous*/
	//InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_water(water_wall_contact, water_air_complex.getInnerRelation());
	//InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration_air(air_wall_contact, air_water_complex.getInnerRelation());
	///** Pressure relaxation for water by using position verlet time stepping. */
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
	//	pressure_relaxation_water(water_wall_contact, water_air_complex);
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
	//	density_relaxation_water(water_wall_contact, water_air_complex);
	///** Extend Pressure relaxation is used for air. */
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
	//	pressure_relaxation_air(air_wall_contact, air_water_complex);
	///*Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
	//	pressure_relaxation_air(air_wall_contact, air_water_complex, 2.0);*/
	//Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
	//	density_relaxation_air(air_wall_contact, air_water_complex);
	//water_block.addBodyStateForRecording<Real>("Pressure");		   // output for debug
	//water_block.addBodyStateForRecording<int>("SurfaceIndicator"); // output for debug
	//air_block.addBodyStateForRecording<Real>("Pressure");		   // output for debug
	//air_block.addBodyStateForRecording<int>("SurfaceIndicator"); // output for debug

	/***Emitter and Disposer ***/
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	/*below is orginal case*/
	/*Vec2d emitter_1_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_1_translation = Vec2d(-DL_sponge-BW, 0.5 * DH) + emitter_1_halfsize;
	BodyAlignedBoxByParticle emitter_1(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(emitter_1_translation)), emitter_1_halfsize));
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_1_inflow_injection(emitter_1, 10, 0);

	Vec2d emitter_2_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_2_translation = Vec2d(-DL_sponge - BW, DH + DH_i) + emitter_2_halfsize;
	BodyAlignedBoxByParticle emitter_2(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(emitter_2_translation)), emitter_2_halfsize));
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_2_inflow_injection(emitter_2, 10, 0);


	Vec2d inlet_buffer_1_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d inlet_buffer_1_translation = Vec2d(-DL_sponge - BW, 0.5 * DH) + inlet_buffer_1_halfsize;
	BodyAlignedBoxByCell emitter_buffer_1(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(inlet_buffer_1_translation)), inlet_buffer_1_halfsize));
	SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity1>> emitter_buffer_1_inflow_condition(emitter_buffer_1);

	Vec2d inlet_buffer_2_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d inlet_buffer_2_translation = Vec2d(-DL_sponge - BW, DH + DH_i) + inlet_buffer_2_halfsize;
	BodyAlignedBoxByCell emitter_buffer_2(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(inlet_buffer_2_translation)), inlet_buffer_2_halfsize));
	SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity1>> emitter_buffer_2_inflow_condition(emitter_buffer_2);
	*/

	/*Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d disposer_translation = Vec2d(DL1 -BW, 0.0) + disposer_halfsize;
	BodyAlignedBoxByCell disposer_water(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(disposer_translation)), disposer_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_water(disposer_water, xAxis);*/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*below is two injection from left*/
	/*for single phase flow*/
	Vec2d emitter_2_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_2_translation = Vec2d(-DL_sponge - BW, DH + DH_i) + emitter_2_halfsize;
	BodyAlignedBoxByParticle emitter_water_2(
		water_block, makeShared<AlignedBoxShape>(Transform2d(emitter_2_translation), emitter_2_halfsize));
	SimpleDynamics<InletInflowCondition> inflow_condition_water_2(emitter_water_2);
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_water_2(emitter_water_2, 350, 0);

	Vec2d emitter_1_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_1_translation = Vec2d(-DL_sponge - BW, 0.5 * DH) + emitter_1_halfsize;
	BodyAlignedBoxByParticle emitter_water_1(
		water_block, makeShared<AlignedBoxShape>(Transform2d(emitter_1_translation), emitter_1_halfsize));
	SimpleDynamics<InletInflowCondition> inflow_condition_water_1(emitter_water_1);
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_water_1(emitter_water_1, 350, 0);

	Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.5 * (1.5 * DH + 2 * DH_i));
	Vec2d disposer_translation = Vec2d(DL1 - BW, 0.0) + disposer_halfsize;
	BodyAlignedBoxByCell disposer_water(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(disposer_translation)), disposer_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_water(disposer_water, xAxis);

	/*for two phases flow*/
	/*Vec2d emitter_2_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_2_translation = Vec2d(-DL_sponge - BW, DH + DH_i) + emitter_2_halfsize;
	BodyAlignedBoxByParticle emitter_water_2(
		water_block, makeShared<AlignedBoxShape>(Transform2d(emitter_2_translation), emitter_2_halfsize));
	SimpleDynamics<InletInflowCondition> inflow_condition_water_2(emitter_water_2);
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_water_2(emitter_water_2, 350, 0);

	Vec2d emitter_1_halfsize = Vec2d(0.5 * BW, 0.5 * DH_i);
	Vec2d emitter_1_translation = Vec2d(-DL_sponge - BW, 0.5 * DH) + emitter_1_halfsize;
	BodyAlignedBoxByParticle emitter_air_1(
		air_block, makeShared<AlignedBoxShape>(Transform2d(emitter_1_translation), emitter_1_halfsize));
	SimpleDynamics<InletInflowCondition> inflow_condition_air_1(emitter_air_1);
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_injection_air_1(emitter_air_1, 350, 0);

	Vec2d disposer_halfsize_water = Vec2d(0.5 * BW, 0.5 * (1.5 * DH + 2 * DH_i));
	Vec2d disposer_translation_water = Vec2d(DL1 - BW, 0.0) + disposer_halfsize_water;
	BodyAlignedBoxByCell disposer_water(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(disposer_translation_water)), disposer_halfsize_water));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_water(disposer_water, xAxis);

	Vec2d disposer_halfsize_air = Vec2d(0.5 * BW, 0.5 * (1.5 * DH + 2 * DH_i));
	Vec2d disposer_translation_air = Vec2d(DL1 - BW, 0.0) + disposer_halfsize_air;
	BodyAlignedBoxByCell disposer_air(
		air_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(disposer_translation_air)), disposer_halfsize_air));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_air(disposer_air, xAxis);*/
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	Real end_time = 100.0;
	Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;							 /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_body_states.writeToFile();
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step_water.exec();
			//initialize_a_fluid_step_air.exec();
			/*for single phase flow*/
			Real Dt = get_fluid_advection_time_step_size_water.exec();
			/*for two phase flow*/
			/*Real Dt_w = get_fluid_advection_time_step_size_water.exec();
			Real Dt_a = get_fluid_advection_time_step_size_air.exec();
			Real Dt = SMIN(Dt_w, Dt_a);*/

			inlet_outlet_surface_particle_indicator_water.exec();
			//inlet_outlet_surface_particle_indicator_air.exec();
			update_density_by_summation_water.exec();
			//update_density_by_summation_air.exec();
			viscous_acceleration_water.exec();
			//viscous_acceleration_air.exec();
			transport_velocity_correction_water.exec();
			//transport_velocity_correction_air.exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/*for single phase flow*/
				dt = SMIN(get_fluid_time_step_size_water.exec(), Dt - relaxation_time);
				/*for two phase flow*/
				/*Real dt_w = get_fluid_time_step_size_water.exec();
				Real dt_a = get_fluid_time_step_size_air.exec();
				dt = SMIN(SMIN(dt_w, dt_a), Dt);*/

				pressure_relaxation_water.exec(dt);
				//pressure_relaxation_air.exec(dt);

				/*main flow with injection on one side*/
				/*for single phase flow*/
				//emitter_buffer_1_inflow_condition.exec();
				//emitter_buffer_2_inflow_condition.exec();
				

				/*two injections on left*/
				/*for single phase flow*/
				inflow_condition_water_1.exec();
				inflow_condition_water_2.exec();
				/*for two phase flow*/
				/*inflow_condition_air_1.exec();
				inflow_condition_water_2.exec();*/

				density_relaxation_water.exec(dt);
				//density_relaxation_air.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			/** inflow injection*/
			/*main flow with injection on one side*/
			/*for single phase flow*/
			//emitter_1_inflow_injection.exec();
			//emitter_2_inflow_injection.exec();
			//disposer_down_outflow_deletion.exec();

			/*two injections on left*/
			/*for single phase flow*/
			emitter_injection_water_1.exec();
			emitter_injection_water_2.exec();
			disposer_outflow_deletion_water.exec();
			/*for two phases flow*/
			/*emitter_injection_air_1.exec();
			emitter_injection_water_2.exec();
			disposer_outflow_deletion_water.exec();
			disposer_outflow_deletion_air.exec();*/

			/** Update cell linked list and configuration. */
			/*for single phase flow*/
			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_complex_relation.updateConfiguration();

			/*for two phases flow*/
			/*water_block.updateCellLinkedListWithParticleSort(100);
			water_air_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();
			air_block.updateCellLinkedListWithParticleSort(100);
			air_water_complex.updateConfiguration();
			air_wall_contact.updateConfiguration();*/
		}

		TickCount t2 = TickCount::now();
		write_body_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;

	return 0;
}
