/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the first geometry.
//	To use this, please commenting the setting for the second geometry.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/James-body.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters for human body.
//----------------------------------------------------------------------
//Vec3d domain_lower_bound(-900.0, -1250.0, -550.0);
//Vec3d domain_upper_bound(900.0, 950.0, 200.0);
//Vecd translation_body(0.0, 0.0, 0.0);
//Real scaling_body = 1.0; 
//----------------------------------------------------------------------
//	Basic geometry parameters for bed.
//----------------------------------------------------------------------
//Real DL1 = 900.0; // Bed length upside.
//Real DL2 = 1200.0; // Bed length downside
//Real DH = 300.0;	// Bed height.
//Real DW = 800.0;	// Bed width half.
//Vecd translation_bed(0.5 * (DL2 - DL1), 0.0, -200.0 - 0.5 * DH);
//Real scaling_bed = 1.0;

/*below for test*/
//Real LL = 1600.0;				  // width
//Real LH = 2300.0;				  // length
//Real LW = 500.0;                  // heigh

/*scaled*/
Vec3d domain_lower_bound(-0.9, -1.25, -0.55);
Vec3d domain_upper_bound(0.9, 0.95, 0.2);
Vecd translation_body(0.0, 0.0, 0.0);
Real scaling_body = 0.001;

/*below for test*/
Real LL = 1.6;				  // width
Real LH = 2.3;				  // length
Real LW = 0.5;                  // heigh
Real scaling_bed = 1.0;

//----------------------------------------------------------------------
// For material properties of the solid.
//----------------------------------------------------------------------
Real gravity_g = 0.9;
Real rho0_s = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5.0e6;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
//----------------------------------------------------------------------
//	define the human body model.
//----------------------------------------------------------------------
class HumanBody : public ComplexShape
{
public:
	explicit HumanBody(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_file, translation_body, scaling_body);
	}
};
//----------------------------------------------------------------------
//	define the bed model.
//----------------------------------------------------------------------
class Bed : public ComplexShape
{
public:
	//explicit Bed (const std::string& shape_name) : ComplexShape(shape_name)
	//{
	//	//Vecd halfsize_bed(0.5 * (DL2 + DL1), DW, 150);
	//	Vecd halfsize_bed();
	//	Transformd translation_bed(translation_bed);
	//	add<TransformShape<GeometricShapeBox>>(Transformd(translation_bed), halfsize_bed);
	//}

	explicit Bed(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_bed(0.5 * LL, 0.5 * LH, 0.5 * LW);
		Transformd translation_bed(Vecd(0.0, -0.1, -0.45));
		add<TransformShape<GeometricShapeBox>>(Transformd(translation_bed), halfsize_bed);
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	IOEnvironment io_environment(system);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	system.setRunParticleRelaxation(false);
	/** Tag for starting with relaxed body-fitted particles distribution */
	system.setReloadParticles(true);

	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	RealBody human_body(system, makeShared<HumanBody>("HumanBody"));
	// level set shape is used for particle relaxation
	human_body.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet()->writeLevelSet(io_environment);
	human_body.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? human_body.generateParticles<ParticleGeneratorReload>(io_environment, human_body.getName())
		: human_body.generateParticles<ParticleGeneratorLattice>();
	human_body.addBodyStateForRecording<Vecd>("PriorAcceleration");
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map.
		//----------------------------------------------------------------------
		InnerRelation imported_model_inner(human_body);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(human_body);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
		//----------------------------------------------------------------------
		//	Define simple file input and outputs functions.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_imported_model_to_vtp(io_environment, { human_body });
		MeshRecordingToPlt write_cell_linked_list(io_environment, human_body.getCellLinkedList());
		ReloadParticleIO write_particle_reload_files(io_environment, { &human_body});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_imported_model_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		write_imported_model_to_vtp.writeToFile(0.0);
		human_body.updateCellLinkedList();
		write_cell_linked_list.writeToFile(0.0);
		//----------------------------------------------------------------------
		//	Particle relaxation time stepping start here.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
				write_imported_model_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of imported model finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Creating bed body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody bed(system, makeShared<Bed>("Bed"));
	bed.defineParticlesAndMaterial<SolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	bed.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, system.real_bodies_);
	//body_states_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation human_body_inner(human_body);
	SurfaceContactRelation human_body_contact(human_body, { &bed });
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> human_body_initialize_timestep(human_body, gravity_ptr);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> human_body_corrected_configuration(human_body_inner);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> human_body_get_time_step_size(human_body);
	/** stress relaxation for the human body. */
	Dynamics1Level<solid_dynamics::Integration1stHalf> human_body_stress_relaxation_first_half(human_body_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> human_body_stress_relaxation_second_half(human_body_inner);
	/** Algorithms for solid-solid contact. */
	InteractionDynamics<solid_dynamics::ContactDensitySummation> human_body_update_contact_density(human_body_contact);
	InteractionDynamics<solid_dynamics::ContactForceFromWall> human_body_compute_solid_contact_forces(human_body_contact);
	/** Damping with the solid body*/
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
		human_body_damping(0.1, human_body_inner, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	human_body_corrected_configuration.exec();
	//----------------------------------------------------------------------
	//	Initial states output.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile(0);
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
	TickCount t1 = TickCount::now();
	TimeInterval interval;
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
				human_body_initialize_timestep.exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				human_body_update_contact_density.exec();
				human_body_compute_solid_contact_forces.exec();
				human_body_stress_relaxation_first_half.exec(dt);
				human_body_damping.exec();
				human_body_stress_relaxation_second_half.exec(dt);

				human_body.updateCellLinkedList();
				human_body_contact.updateConfiguration();

				ite++;
				dt = human_body_get_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}
		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
