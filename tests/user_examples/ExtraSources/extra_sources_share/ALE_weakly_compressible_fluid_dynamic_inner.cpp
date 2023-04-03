#include "ALE_weakly_compressible_fluid_dynamic_inner.h"
#include "ALE_weakly_compressible_fluid_dynamic_inner.hpp"


namespace SPH
{
	//=================================================================================================//
	namespace ALE_weakly_compressible_fluid_dynamic
	{
		ALEFlowTimeStepInitialization::
			ALEFlowTimeStepInitialization(SPHBody& sph_body, SharedPtr<Gravity> gravity_ptr)
			: BaseTimeStepInitialization(sph_body, gravity_ptr),
			EulerianWeaklyCompressibleFluidDataSimple(sph_body), rho_(particles_->rho_),
			pos_(particles_->pos_), dmom_dt_prior_(particles_->dmom_dt_prior_) {}
		//=================================================================================================//
		void ALEFlowTimeStepInitialization::update(size_t index_i, Real dt)
		{
			dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
		}
		//=================================================================================================//
		BaseDensitySummationInner::BaseDensitySummationInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			rho_(particles_->rho_), rho_sum_(particles_->rho_sum_), mass_(particles_->mass_),
			rho0_(sph_body_.base_material_->ReferenceDensity()) {}
		//=================================================================================================//
		void BaseDensitySummationInner::update(size_t index_i, Real dt)
		{
			rho_[index_i] = rho_sum_[index_i];
		}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseInnerRelation& inner_relation)
			: BaseDensitySummationInner(inner_relation),
			W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)),
			inv_sigma0_(1.0 / sph_body_.sph_adaptation_->ReferenceNumberDensity()) {}
		//=================================================================================================//
		void DensitySummationInner::interaction(size_t index_i, Real dt)
		{
			Real sigma = W0_;
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
		}
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseInnerRelation& inner_relation, Real coefficient)
			: LocalDynamics(inner_relation.getSPHBody()), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			coefficient_(coefficient) {}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration_trans = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

				// acceleration for transport velocity
				acceleration_trans -= 2.0 * nablaW_ijV_j;
			}

			if (surface_indicator_[index_i] == 0)
				pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
		}
		//=================================================================================================//
		TransportVelocityCorrectionInnerAdaptive::
			TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation& inner_relation, Real coefficient)
			: LocalDynamics(inner_relation.getSPHBody()), EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			sph_adaptation_(*sph_body_.sph_adaptation_),
			pos_(particles_->pos_), surface_indicator_(particles_->surface_indicator_),
			smoothing_length_sqr_(pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
			coefficient_(coefficient) {}
		//=================================================================================================//
		void TransportVelocityCorrectionInnerAdaptive::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration_trans = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

				// acceleration for transport velocity
				acceleration_trans -= 2.0 * nablaW_ijV_j;
			}

			if (surface_indicator_[index_i] == 0)
			{
				Real inv_h_ratio = 1.0 / sph_adaptation_.SmoothingLengthRatio(index_i);
				pos_[index_i] += coefficient_ * smoothing_length_sqr_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
			}
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation),
			Vol_(particles_->Vol_), rho_(particles_->rho_), p_(particles_->p_),
			vel_(particles_->vel_), dmom_dt_prior_(particles_->dmom_dt_prior_),
			mu_(particles_->fluid_.ReferenceViscosity()),
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			const Vecd& vel_i = vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				// viscous force
				vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
			}

			dmom_dt_prior_[index_i] += rho_i * acceleration;
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(SPHBody& sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
			EulerianWeaklyCompressibleFluidDataSimple(sph_body),
			fluid_(particles_->fluid_), rho_(particles_->rho_),
			p_(particles_->p_), vel_(particles_->vel_),
			smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::outputResult(Real reduced_value)
		{
			// since the particle does not change its configuration in pressure relaxation step
			// I chose a time-step size according to Eulerian method
			return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSizeForImplicitViscosity::
			AdvectionTimeStepSizeForImplicitViscosity(SPHBody& sph_body, Real U_max, Real advectionCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max* U_max),
			EulerianWeaklyCompressibleFluidDataSimple(sph_body), vel_(particles_->vel_),
			smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			advectionCFL_(advectionCFL) {}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::reduce(size_t index_i, Real dt)
		{
			return vel_[index_i].squaredNorm();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSizeForImplicitViscosity::outputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			return advectionCFL_ * smoothing_length_min_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL)
			: AdvectionTimeStepSizeForImplicitViscosity(sph_body, U_max, advectionCFL),
			fluid_(particles_->fluid_)
		{
			Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
			reference_ = SMAX(viscous_speed * viscous_speed, reference_);
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return AdvectionTimeStepSizeForImplicitViscosity::reduce(index_i, dt);
		}
		//=================================================================================================//
		BaseIntegration::BaseIntegration(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()),
			EulerianWeaklyCompressibleFluidDataInner(inner_relation), fluid_(particles_->fluid_),
			Vol_(particles_->Vol_), mass_(particles_->mass_), rho_(particles_->rho_),
			p_(particles_->p_), drho_dt_(particles_->drho_dt_), vel_(particles_->vel_), mom_(particles_->mom_),
			dmom_dt_(particles_->dmom_dt_), dmom_dt_prior_(particles_->dmom_dt_prior_), pos_(particles_->pos_) {}
		//=================================================================================================//
		FluidStarState AcousticRiemannSolverOLD::
			getInterfaceState(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
		{
			Real ul = -e_ij.dot(state_i.vel_);
			Real ur = -e_ij.dot(state_j.vel_);
			Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
			Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
			Real p_star = 0;
			Vecd v_star = Vecd::Zero();

			Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);
			p_star = (rhol_cl * state_j.p_ + rhor_cr * state_j.p_ + rhol_cl * rhor_cr * (ul - ur) * SMIN(3.0 * SMAX((ul - ur) / clr, 0.0), 1.0)) /
				(rhol_cl + rhor_cr);
			Real u_star = (rhol_cl * ul + rhor_cr * ur + state_i.p_ - state_j.p_) / (rhol_cl + rhor_cr);
			v_star = (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) -
				e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));

			FluidStarState interface_state(v_star, p_star);

			return interface_state;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//