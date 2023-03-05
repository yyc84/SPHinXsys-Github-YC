/**
 * @file 	structural_simulation_class.h
 * @brief 	The structural simulation module is licensed under the Aladdin Free Public License (https://spdx.org/licenses/Aladdin.html) regarding usage for medical device development.
 * Commercial use for medical device development is not permitted. This does not apply to applications in other fields.
 * @details	solid structural simulation class for general structural simulations
 * @author 	Bence Z. Rochlitz - Virtonomy GmbH, Xiangyu Hu
 */

#ifndef SOLID_STRUCTURAL_SIMULATION_CLASS_H
#define SOLID_STRUCTURAL_SIMULATION_CLASS_H

#include "sphinxsys.h"
#include <algorithm>
#include <memory>
#include <vector>

using namespace SPH;
using namespace std;
using GravityPair = pair<int, Vec3d>;
using AccelTuple = tuple<int, BoundingBox, Vec3d>;
using ForceTuple = tuple<int, BoundingBox, Vec3d, Real>;
using PressureTuple = tuple<int, SharedPtr<TriangleMeshShape>, Vec3d, StdVec<array<Real, 2>>>;
using SpringDamperTuple = tuple<int, Vec3d, Real>;
/**
 * @brief SurfaceSpringTuple
 * int: body index
 * TriangleMeshShape*: the body part, the normal spring is applied to
 * bool: if true, the "outer" surface is considered (particle normals > 90° from the particle-source point vector), if false, the "inner" surface
 * Vec3d: source point to relate inner and outer surface
 * Real: normal spring stiffness
 * Real: damping coefficient
 */
using SurfaceSpringTuple = tuple<int, SharedPtr<TriangleMeshShape>, bool, Vec3d, Real, Real>;
using ConstrainedRegionPair = pair<int, BoundingBox>;
using PositionSolidBodyTuple = tuple<int, Real, Real, Vec3d>;
using PositionScaleSolidBodyTuple = tuple<int, Real, Real, Real>;
using TranslateSolidBodyTuple = tuple<int, Real, Real, Vec3d>;
using TranslateSolidBodyPartTuple = tuple<int, Real, Real, Vec3d, BoundingBox>;

#ifdef __EMSCRIPTEN__
struct StlData
{
	string name;
	uintptr_t ptr;
};

using StlList = vector<StlData>;
#else
using StlList = vector<string>;
#endif

class BodyPartFromMesh : public BodyRegionByParticle
{
public:
	BodyPartFromMesh(SPHBody &body, SharedPtr<TriangleMeshShape> triangle_mesh_shape_ptr);
	~BodyPartFromMesh(){};
};

class SolidBodyFromMesh : public SolidBody
{
public:
	SolidBodyFromMesh(SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
					  SharedPtr<SaintVenantKirchhoffSolid> material_model, StdLargeVec<Vec3d> &pos_0, StdLargeVec<Real> &volume);
	~SolidBodyFromMesh(){};
};

class SolidBodyForSimulation
{
private:
	SolidBodyFromMesh solid_body_from_mesh_;
	InnerRelation inner_body_relation_;

	SimpleDynamics<NormalDirectionFromBodyShape> initial_normal_direction_;
	InteractionDynamics<solid_dynamics::CorrectConfiguration> correct_configuration_;
	Dynamics1Level<solid_dynamics::Integration1stHalf> stress_relaxation_first_half_;
	Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>> damping_random_;

public:
	// no particle reload --> direct generator
	SolidBodyForSimulation(
		SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
		Real physical_viscosity, SharedPtr<SaintVenantKirchhoffSolid> material_model, StdLargeVec<Vec3d> &pos_0, StdLargeVec<Real> &volume);
	~SolidBodyForSimulation(){};

	SolidBodyFromMesh *getSolidBodyFromMesh() { return &solid_body_from_mesh_; };
	ElasticSolidParticles *getElasticSolidParticles() { return DynamicCast<ElasticSolidParticles>(this, &solid_body_from_mesh_.getBaseParticles()); };
	InnerRelation *getInnerBodyRelation() { return &inner_body_relation_; };

	SimpleDynamics<NormalDirectionFromBodyShape> *getInitialNormalDirection() { return &initial_normal_direction_; };
	InteractionDynamics<solid_dynamics::CorrectConfiguration> *getCorrectConfiguration() { return &correct_configuration_; };
	Dynamics1Level<solid_dynamics::Integration1stHalf> *getStressRelaxationFirstHalf() { return &stress_relaxation_first_half_; };
	Dynamics1Level<solid_dynamics::Integration2ndHalf> *getStressRelaxationSecondHalf() { return &stress_relaxation_second_half_; };
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>> *getDampingWithRandomChoice() { return &damping_random_; };
};

void expandBoundingBox(BoundingBox *original, BoundingBox *additional);

void relaxParticlesSingleResolution(IOEnvironment &io_environment,
									bool write_particles_to_file,
									SolidBodyFromMesh &solid_body_from_mesh,
									ElasticSolidParticles &solid_body_from_mesh_particles,
									InnerRelation &solid_body_from_mesh_inner);

static inline Real getPhysicalViscosityGeneral(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 1.0)
{
	// the physical viscosity is defined in the paper pf prof. Hu
	// https://arxiv.org/pdf/2103.08932.pdf
	// physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
	// beta: shape constant --> how to define it? - it's 1 for now.. TODO
	// L: length scale of the problem --> 10 mm roughly
	return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

class StructuralSimulationInput
{
public:
	string relative_input_path_;
	StlList imported_stl_list_;
	Real scale_stl_;
	vector<Vec3d> translation_list_;
	vector<Real> resolution_list_;
	vector<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list_;
	StdVec<Real> physical_viscosity_;
	StdVec<IndexVector> contacting_body_pairs_list_;
	vector<pair<array<int, 2>, array<Real, 2>>> time_dep_contacting_body_pairs_list_;
	// scale system boundaries
	Real scale_system_boundaries_;
	// particle relaxation
	vector<bool> particle_relaxation_list_;
	bool write_particle_relaxation_data_;
	// boundary conditions
	vector<GravityPair> non_zero_gravity_;
	vector<AccelTuple> acceleration_bounding_box_tuple_;
	vector<ForceTuple> force_in_body_region_tuple_;
	vector<PressureTuple> surface_pressure_tuple_;
	vector<SpringDamperTuple> spring_damper_tuple_;
	vector<SurfaceSpringTuple> surface_spring_tuple_;
	vector<int> body_indices_fixed_constraint_;
	vector<ConstrainedRegionPair> body_indices_fixed_constraint_region_;
	vector<PositionSolidBodyTuple> position_solid_body_tuple_;
	vector<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
	vector<TranslateSolidBodyTuple> translation_solid_body_tuple_;
	vector<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

	StructuralSimulationInput(
		string relative_input_path,
		StlList imported_stl_list,
		Real scale_stl,
		vector<Vec3d> translation_list,
		vector<Real> resolution_list,
		vector<shared_ptr<SaintVenantKirchhoffSolid>> material_model_list,
		StdVec<Real> physical_viscosity,
		StdVec<IndexVector> contacting_bodies_list);
};

class StructuralSimulation
{
private:
	UniquePtrKeepers<SurfaceContactRelation> contact_relation_ptr_keeper_;
	UniquePtrKeepers<Gravity> gravity_ptr_keeper_;
	UniquePtrKeepers<BodyPartFromMesh> body_part_tri_mesh_ptr_keeper_;

protected:
	// mandatory input
	string relative_input_path_;
	StlList imported_stl_list_;
	Real scale_stl_;
	vector<Vec3d> translation_list_;
	vector<Real> resolution_list_;
	vector<SharedPtr<TriangleMeshShape>> body_mesh_list_;
	vector<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list_;
	StdVec<Real> physical_viscosity_;
	StdVec<IndexVector> contacting_body_pairs_list_;
	vector<pair<array<int, 2>, array<Real, 2>>> time_dep_contacting_body_pairs_list_; // optional: time dependent contact
	vector<bool> particle_relaxation_list_;											  // optional: particle relaxation
	bool write_particle_relaxation_data_;

	// internal members
	Real system_resolution_;
	SPHSystem system_;
	Real scale_system_boundaries_;
	IOEnvironment io_environment_;

	vector<shared_ptr<SolidBodyForSimulation>> solid_body_list_;
	vector<shared_ptr<SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection>>> particle_normal_update_;

	vector<shared_ptr<SurfaceContactRelation>> contact_list_;
	vector<shared_ptr<InteractionDynamics<solid_dynamics::ContactDensitySummation>>> contact_density_list_;
	vector<shared_ptr<InteractionDynamics<solid_dynamics::ContactForce>>> contact_force_list_;

	// for initializeATimeStep
	vector<shared_ptr<SimpleDynamics<TimeStepInitialization>>> initialize_time_step_;
	vector<GravityPair> non_zero_gravity_;
	// for AccelerationForBodyPartInBoundingBox
	vector<shared_ptr<SimpleDynamics<solid_dynamics::AccelerationForBodyPartInBoundingBox>>> acceleration_bounding_box_;
	vector<AccelTuple> acceleration_bounding_box_tuple_;
	// for ForceInBodyRegion
	vector<shared_ptr<SimpleDynamics<solid_dynamics::ForceInBodyRegion>>> force_in_body_region_;
	vector<ForceTuple> force_in_body_region_tuple_;
	// for SurfacePressureFromSource
	vector<shared_ptr<SimpleDynamics<solid_dynamics::SurfacePressureFromSource>>> surface_pressure_;
	vector<PressureTuple> surface_pressure_tuple_;
	// for SpringDamperConstraintParticleWise
	vector<shared_ptr<SimpleDynamics<solid_dynamics::SpringDamperConstraintParticleWise>>> spring_damper_constraint_;
	vector<SpringDamperTuple> spring_damper_tuple_;
	// for SpringNormalOnSurfaceParticles
	vector<shared_ptr<SimpleDynamics<solid_dynamics::SpringNormalOnSurfaceParticles>>> surface_spring_;
	vector<SurfaceSpringTuple> surface_spring_tuple_;
	// for ConstrainSolidBody
	vector<shared_ptr<SimpleDynamics<solid_dynamics::FixBodyConstraint>>> fixed_constraint_body_;
	vector<int> body_indices_fixed_constraint_;
	// for ConstrainSolidBodyRegion
	vector<shared_ptr<SimpleDynamics<solid_dynamics::FixBodyPartConstraint>>> fixed_constraint_region_;
	vector<ConstrainedRegionPair> body_indices_fixed_constraint_region_;
	// for PositionSolidBody
	vector<shared_ptr<SimpleDynamics<solid_dynamics::PositionSolidBody>>> position_solid_body_;
	vector<PositionSolidBodyTuple> position_solid_body_tuple_;
	// for PositionScaleSolidBody
	vector<shared_ptr<SimpleDynamics<solid_dynamics::PositionScaleSolidBody>>> position_scale_solid_body_;
	vector<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
	// for TranslateSolidBody
	vector<shared_ptr<SimpleDynamics<solid_dynamics::TranslateSolidBody>>> translation_solid_body_;
	vector<TranslateSolidBodyTuple> translation_solid_body_tuple_;
	// for TranslateSolidBodyPart
	vector<shared_ptr<SimpleDynamics<solid_dynamics::TranslateSolidBodyPart>>> translation_solid_body_part_;
	vector<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

	// iterators
	int iteration_;

	// data storage
	vector<Real> von_mises_stress_max_;
	StdLargeVec<StdLargeVec<Real>> von_mises_stress_particles_;

	vector<Real> von_mises_strain_max_;
	StdLargeVec<StdLargeVec<Real>> von_mises_strain_particles_;

	// for constructor, the order is important
	void scaleTranslationAndResolution();
	void setSystemResolutionMax();
	void createBodyMeshList();
	void calculateSystemBoundaries();
	void initializeElasticSolidBodies();
	void initializeContactBetweenTwoBodies(int first, int second);
	void initializeAllContacts();

	// for initializeBoundaryConditions
	void initializeGravity();
	void initializeAccelerationForBodyPartInBoundingBox();
	void initializeForceInBodyRegion();
	void initializeSurfacePressure();
	void initializeSpringDamperConstraintParticleWise();
	void initializeSpringNormalOnSurfaceParticles();
	void initializeConstrainSolidBody();
	void initializeConstrainSolidBodyRegion();
	void initializePositionSolidBody();
	void initializePositionScaleSolidBody();
	void initializeTranslateSolidBody();
	void initializeTranslateSolidBodyPart();

	// for runSimulation, the order is important
	void executeInitialNormalDirection();
	void executeCorrectConfiguration();
	void executeUpdateElasticNormalDirection();
	void executeInitializeATimeStep();
	void executeAccelerationForBodyPartInBoundingBox();
	void executeForceInBodyRegion();
	void executeSurfacePressure();
	void executeSpringDamperConstraintParticleWise();
	void executeSpringNormalOnSurfaceParticles();
	void executeContactDensitySummation();
	void executeContactForce();
	void executeStressRelaxationFirstHalf(Real dt);
	void executeConstrainSolidBody();
	void executeConstrainSolidBodyRegion();
	void executePositionSolidBody(Real dt);
	void executePositionScaleSolidBody(Real dt);
	void executeTranslateSolidBody(Real dt);
	void executeTranslateSolidBodyPart(Real dt);
	void executeDamping(Real dt);
	void executeStressRelaxationSecondHalf(Real dt);
	void executeUpdateCellLinkedList();
	void executeContactUpdateConfiguration();

	void initializeSimulation();

	void runSimulationStep(Real &dt, Real &integration_time);

public:
	explicit StructuralSimulation(const StructuralSimulationInput &input);
	~StructuralSimulation();

	StdVec<shared_ptr<SolidBodyForSimulation>> get_solid_body_list_() { return solid_body_list_; };
	Real getMaxDisplacement(int body_index);

	// For c++
	void runSimulation(Real end_time);

	// For JS
	double runSimulationFixedDurationJS(int number_of_steps);
};

class StructuralSimulationJS : public StructuralSimulation
{
public:
	StructuralSimulationJS(const StructuralSimulationInput &input);
	~StructuralSimulationJS() = default;

	void runSimulationFixedDuration(int number_of_steps);

	VtuStringData getVtuData();

private:
	BodyStatesRecordingToVtpString write_states_;
	Real dt;
};

#endif // SOLID_STRUCTURAL_SIMULATION_CLASS_H
