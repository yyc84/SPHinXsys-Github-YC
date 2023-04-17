/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	particle_dynamics_diffusion_reaction.h
 * @brief 	This is the particle dynamics applicable for all type bodies
 * 			TODO: there is an issue on applying corrected configuration for contact bodies.. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef EXTRA_TWO_PHASE_HEAT_TRANSFER_H
#define EXTRA_TWO_PHASE_HEAT_TRANSFER_H

#include "all_particle_dynamics.h"
#include "two_phase_heat_transfer_particles.h"
#include "diffusion_reaction.h"

namespace SPH
{
	template <class DiffusionReactionParticlesType>
	using TwoPhaseDiffusionReactionSimpleData =
		DataDelegateSimple<DiffusionReactionParticlesType>;

	template <class DiffusionReactionParticlesType>
	using TwoPhaseDiffusionReactionInnerData =
		DataDelegateInner< DiffusionReactionParticlesType>;

	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	using TwoPhaseDiffusionReactionContactData =
		DataDelegateContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType,
		DataDelegateEmptyBase>;
	/**
	 * @class  DiffusionReactionInitialCondition
	 * @brief pure abstract class for initial conditions
	 */
	template <class DiffusionReactionParticlesType>
	class TwoPhaseDiffusionReactionInitialCondition
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	public:
		explicit TwoPhaseDiffusionReactionInitialCondition(SPHBody &sph_body);
		virtual ~TwoPhaseDiffusionReactionInitialCondition(){};

	protected:
		StdLargeVec<Vecd> &pos_;
		StdVec<StdLargeVec<Real>>& all_species_;
		StdLargeVec<Real>& thermal_conductivity_;
	};

	/**
	 * @class GetDiffusionTimeStepSize
	 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
	 */
	template <class DiffusionReactionParticlesType>
	class TwoPhaseGetDiffusionTimeStepSize
		: public BaseDynamics<Real>,
		  public TwoPhaseDiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	public:
		explicit TwoPhaseGetDiffusionTimeStepSize(SPHBody &sph_body);
		virtual ~TwoPhaseGetDiffusionTimeStepSize(){};

		virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };

	protected:
		Real diff_time_step_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesInner
	 * @brief Compute the diffusion relaxation process of all species
	 */
	template <class DiffusionReactionParticlesType>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesInner
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionInnerData<DiffusionReactionParticlesType>
	{
	protected:
		typedef typename DiffusionReactionParticlesType::DiffusionReactionMaterial DiffusionReactionMaterial;
		DiffusionReactionMaterial & material_;
		StdVec<BaseDiffusion*> &all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>*>& gradient_species_;
		StdVec<StdLargeVec<Real>*> diffusion_dt_;

		void initializeDiffusionChangeRate(size_t particle_i);
		void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij);
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

	public:
		typedef DiffusionReactionParticlesType InnerParticlesType;
		typedef BaseInnerRelation BodyRelationType;

		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesInner(BaseInnerRelation& inner_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesInner() {};
		StdVec<BaseDiffusion*>& AllDiffusions() { return material_.AllDiffusions(); };
		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	};


	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesComplex
		: public TwoPhaseRelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>,
		public TwoPhaseDiffusionReactionContactData<DiffusionReactionParticlesType,
		ContactDiffusionReactionParticlesType>
	{
		StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;
		StdVec<StdVec<StdLargeVec<Real>*>> contact_thermal_conductivity_;
		StdLargeVec<Real> &thermal_conductivity_;
		StdLargeVec<Real> &external_diffusion_dt_;
		StdLargeVec<Real> &external_diffusion_dt_sum_;
	protected:
		void initializeExternalDiffusionChangeRate(size_t particle_i);
		void getDiffusionChangeRateContact(size_t particle_i, size_t particle_j,
			Vecd& e_ij, Real surface_area_ij,
			const StdVec<StdLargeVec<Real>*>& gradient_species_k, const StdVec<StdLargeVec<Real>*>& thermal_conductivity_k);

		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		typedef ComplexRelation BodyRelationType;
		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesComplex(ComplexRelation& complex_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesComplex() {};
		void interaction(size_t index_i, Real dt = 0.0);

	};

	/**
	 * @class InitializationRK
	 * @brief initialization of a runge-kutta integration scheme
	 */
	template <class DiffusionReactionParticlesType>
	class TwoPhaseInitializationRK
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	protected:
		typename DiffusionReactionParticlesType::DiffusionReactionMaterial &material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;

	public:
		TwoPhaseInitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s);
		virtual ~TwoPhaseInitializationRK(){};

		void update(size_t index_i, Real dt = 0.0);
	};

	
	/**
	 * @class SecondStageRK2
	 * @brief the second stage of the 2nd-order Runge-Kutta scheme
	 */
	template <class FirstStageType>
	class TwoPhaseSecondStageRK2 : public FirstStageType
	{

	protected:
		StdVec<StdLargeVec<Real>>&diffusion_species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;
	public:
		TwoPhaseSecondStageRK2(typename FirstStageType::BodyRelationType& body_relation,
			StdVec<StdLargeVec<Real>>& species_s);
		virtual ~TwoPhaseSecondStageRK2() {};
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	template <class FirstStageType>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesRK2 : public BaseDynamics<void>
	{
	protected:
		StdVec<StdLargeVec<Real>> diffusion_species_s_;
		SimpleDynamics<InitializationRK<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<SecondStageRK2<FirstStageType>> rk2_2nd_stage_;
		StdVec<BaseDiffusion*> all_diffusions_;

	public:
		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType& body_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesRK2() {};

		virtual void exec(Real dt = 0.0) override;
	};

}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_H