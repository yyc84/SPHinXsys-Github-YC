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
 * @file 	static_confinement.h
 * @brief 	Here, we define the static confinement boundary condition classes for fluid dynamics.
 * @details     This boundary condition is based on Level-set filed.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef FLUID_BOUNDARY_STATIC_COFINEMENT_H
#define FLUID_BOUNDARY_STATIC_COFINEMENT_H


#include "fluid_dynamics_inner.h"
#include "all_particle_dynamics.h"
#include "all_general_dynamics.h"
#include "base_kernel.h"
#include "external_force.h"
#include "riemann_solver.h"
#include "compressible_fluid.h"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"


namespace SPH
{
	namespace ALE_weakly_compressible_fluid_dynamic
	{
		typedef DataDelegateSimple<WeaklyCompressibleFluidParticles> EulerianWeaklyCompressibleFluidDataSimple;
		typedef DataDelegateInner<WeaklyCompressibleFluidParticles> EulerianWeaklyCompressibleFluidDataInner;
		
		class ALEFlowTimeStepInitialization
			: public BaseTimeStepInitialization,
			  public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
			StdLargeVec<Real> &rho_;
			StdLargeVec<Vecd> &pos_, &dmom_dt_prior_;

		public:
			ALEFlowTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
			virtual ~ALEFlowTimeStepInitialization(){};

			void update(size_t index_i, Real dt = 0.0);
		};
		
		class BaseDensitySummationInner : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit BaseDensitySummationInner(BaseInnerRelation& inner_relation);
			virtual ~BaseDensitySummationInner() {};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			Real rho0_;
			StdLargeVec<Real>& rho_, & rho_sum_, & mass_;
		};

		/**
		 * @class DensitySummationInner
		 * @brief  computing density by summation
		 */
		class DensitySummationInner : public BaseDensitySummationInner
		{
		public:
			explicit DensitySummationInner(BaseInnerRelation& inner_relation);
			virtual ~DensitySummationInner() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real W0_, inv_sigma0_;
		};
		class TransportVelocityCorrectionInner : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit TransportVelocityCorrectionInner(BaseInnerRelation& inner_relation, Real coefficient = 0.2);
			virtual ~TransportVelocityCorrectionInner() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd>& pos_;
			StdLargeVec<int>& surface_indicator_;
			Real smoothing_length_sqr_;
			const Real coefficient_;
		};

		/**
		 * @class TransportVelocityCorrectionInner
		 * @brief transport velocity correction
		 */
		class TransportVelocityCorrectionInnerAdaptive : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation& inner_relation, Real coefficient = 0.2);
			virtual ~TransportVelocityCorrectionInnerAdaptive() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			SPHAdaptation& sph_adaptation_;
			StdLargeVec<Vecd>& pos_;
			StdLargeVec<int>& surface_indicator_;
			Real smoothing_length_sqr_;
			const Real coefficient_;
		};
		class ViscousAccelerationInner
			: public LocalDynamics,
			public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit ViscousAccelerationInner(BaseInnerRelation& inner_relation);
			virtual ~ViscousAccelerationInner() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real>& Vol_, & rho_, & p_;
			StdLargeVec<Vecd>& vel_, & dmom_dt_prior_;
			Real mu_;
			Real smoothing_length_;
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>,
			public EulerianWeaklyCompressibleFluidDataSimple
		{
		protected:
			Fluid& fluid_;
			StdLargeVec<Real>& rho_, & p_;
			StdLargeVec<Vecd>& vel_;
			Real smoothing_length_;

		public:
			explicit AcousticTimeStepSize(SPHBody& sph_body);
			virtual ~AcousticTimeStepSize() {};

			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;
		};
		class AdvectionTimeStepSizeForImplicitViscosity
			: public LocalDynamicsReduce<Real, ReduceMax>,
			public EulerianWeaklyCompressibleFluidDataSimple
		{
		public:
			explicit AdvectionTimeStepSizeForImplicitViscosity(
				SPHBody& sph_body, Real U_max, Real advectionCFL = 0.25);
			virtual ~AdvectionTimeStepSizeForImplicitViscosity() {};
			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;

		protected:
			Real smoothing_length_min_;
			StdLargeVec<Vecd>& vel_;
			Real advectionCFL_;
		};

		class AdvectionTimeStepSize : public AdvectionTimeStepSizeForImplicitViscosity
		{
		public:
			explicit AdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL = 0.25);
			virtual ~AdvectionTimeStepSize() {};
			Real reduce(size_t index_i, Real dt = 0.0);

		protected:
			Fluid& fluid_;
		};
		/** define the Acoustic Riemann solver with old form */
		class AcousticRiemannSolverOLD
		{
			Fluid& fluid_i_, & fluid_j_;
		public:
			AcousticRiemannSolverOLD(Fluid& fluid_i, Fluid& fluid_j) : fluid_i_(fluid_i), fluid_j_(fluid_j) {};
			FluidStarState getInterfaceState(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
		};
		/** define the BaseIntergration with pos */
		class BaseIntegration : public LocalDynamics, public EulerianWeaklyCompressibleFluidDataInner
		{
		public:
			explicit BaseIntegration(BaseInnerRelation& inner_relation);
			virtual ~BaseIntegration() {};

		protected:
			Fluid& fluid_;
			StdLargeVec<Real>& Vol_, & mass_, & rho_, & p_, & drho_dt_;
			StdLargeVec<Vecd>& vel_, & mom_, & dmom_dt_, & dmom_dt_prior_, & pos_;
		};

		/**
		 * @class BaseIntegration1stHalf
		 * @brief Template class for pressure relaxation scheme with the Riemann solver
		 * as template variable
		 */
		template <class RiemannSolverType>
		class BaseIntegration1stHalf : public BaseIntegration
		{
		public:
			explicit BaseIntegration1stHalf(BaseInnerRelation& inner_relation);
			virtual ~BaseIntegration1stHalf() {};
			RiemannSolverType riemann_solver_;
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};
		using Integration1stHalf = BaseIntegration1stHalf<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using Integration1stHalfAcousticRiemann = BaseIntegration1stHalf<AcousticRiemannSolver>;

		using Integration1stHalfHLLCRiemann = BaseIntegration1stHalf<HLLCRiemannSolverInWeaklyCompressibleFluid>;
		using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;

		using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;
		using Integration1stHalfAcousticRiemannOLD = BaseIntegration1stHalf<AcousticRiemannSolverOLD>;

		/**
		 * @class Integration2ndHalf
		 * @brief  Template density relaxation scheme with different Riemann solver
		 */
		template <class RiemannSolverType>
		class BaseIntegration2ndHalf : public BaseIntegration
		{
		public:
			explicit BaseIntegration2ndHalf(BaseInnerRelation& inner_relation);
			virtual ~BaseIntegration2ndHalf() {};
			RiemannSolverType riemann_solver_;
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};
		using Integration2ndHalf = BaseIntegration2ndHalf<NoRiemannSolver>;
		using Integration2ndHalfAcousticRiemann = BaseIntegration2ndHalf<AcousticRiemannSolver>;
		using Integration2ndHalfHLLCWithLimiterRiemann = BaseIntegration2ndHalf<HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid>;
		using Integration2ndHalfAcousticRiemannOLD = BaseIntegration2ndHalf<AcousticRiemannSolverOLD>;

	}
}
#endif FLUID_BOUNDARY_STATIC_COFINEMENT_H