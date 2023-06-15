/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                         *
 * ------------------------------------------------------------------------*/
#pragma once

#include "common_compressible_FVM_classes.h"


namespace SPH
{
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration1stHalfInFVM<RiemannSolverType>::BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressibleFVM(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_, limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::initialization(size_t index_i, Real dt)
	{
		E_[index_i] += dE_dt_[index_i] * dt * 0.5;
		rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		Real rho_e = E_[index_i] - 0.5 * mom_[index_i].squaredNorm() / rho_[index_i];
		p_[index_i] = compressible_fluid_.getPressure(rho_[index_i], rho_e);
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
		Vecd momentum_change_rate = dmom_dt_prior_[index_i];
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			if (inner_neighborhood.boundary_type_[n] == 2)
			{
				CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				momentum_change_rate -= 2.0 * dW_ijV_j *
				((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
			}

			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				//non-slip wall boundary
				/*Vecd vel_in_wall = -state_i.vel_;
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;*/

				//rigid wall boundary 
				Vecd vel_in_wall = (state_i.vel_ - e_ij.dot(state_i.vel_)*(e_ij)) + (-e_ij.dot(state_i.vel_)*(e_ij));
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;
				Real E_in_wall=state_i.E_;

				CompressibleFluidState state_j(rho_in_wall, vel_in_wall, p_in_wall, E_in_wall);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				momentum_change_rate -= 2.0 * dW_ijV_j *
				((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
			}

            if (inner_neighborhood.boundary_type_[n] == 10)
            {
				// parabolic velocity inflow
                //Real U_f = 1.0; //characteristic velocity is set as 1
                //Real h = 1.0;   // the height of the inflow domain size
                //Vecd parabolic_velocity_inlet = Vecd::Zero();
                //parabolic_velocity_inlet[0] = 1.5 * U_f * (1.0 - pos_[index_i][1] * pos_[index_i][1] / pow(0.5 * h, 2));
                //Vecd vel_inlet = parabolic_velocity_inlet;
                //Real p_inlet = state_i.p_;
                //Real rho_inlet = state_i.rho_;

				//given value inlet flow
				Real rho0_another = 8.0;					/**< initial density of another. */
				Real u_another = 8.25*sin(3.14159 / 3.0);	/**< initial velocity of another in X axis. */
				Real v_another = -8.25*cos(3.14159 / 3.0);	/**< initial velocity of another in Y axis. */
				Vecd vel_another= Vecd::Zero();
				vel_another[0] = u_another;
				vel_another[1] = v_another;
				Real p_another = 140.2 / 1.2;					/**< initial pressure of another. */
				Real rho_e_another = p_another / (1.4 - 1.0);
				Real E_inlet_another = rho_e_another + 0.5 * rho0_another * vel_another.squaredNorm();

				Real rho_inlet = rho0_another;
				Real p_inlet = p_another;
				Vecd vel_inlet= Vecd::Zero();
				vel_inlet[0] = u_another;
				vel_inlet[1] = v_another;
				Real E_inlet = E_inlet_another;
                CompressibleFluidState state_j(rho_inlet, vel_inlet, p_inlet,E_inlet);
                CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

                momentum_change_rate -= 2.0 * dW_ijV_j *
				((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
            }
                
            if (inner_neighborhood.boundary_type_[n] == 36)
            {
                //Outlet boundary condition
                Vecd vel_outlet = state_i.vel_;
                Real p_outlet = state_i.p_;
                Real rho_outlet = state_i.rho_;
                Real E_outlet = state_i.E_;
                CompressibleFluidState state_j(rho_outlet, vel_outlet, p_outlet, E_outlet);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

                momentum_change_rate -= 2.0 * dW_ijV_j *
				((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
            }

			//Top boundary condition
            if (inner_neighborhood.boundary_type_[n] == 4)
            {
                Real rho0_one = 1.4;                         /**< initial density of one fluid. */
                Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
                Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
                Real p_one = 1.0;                            /**< initial pressure of one fluid. */
                Real rho0_another = 8.0;                     /**< initial density of another. */
                Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
                Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
                Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
                Real gamma_ = 1.4;

                Real run_time = GlobalStaticVariables::physical_time_;
                Real x_1 = 1.0 / 6.0 + run_time * 10.0 / sin(3.14159 / 3.0);
                Real p_top = 0.0;
                Real rho_top = 0.0;
                Vecd vel_top = Vecd::Zero();
                Real E_top = 0.0;
                if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - x_1))
                {
                    rho_top = rho0_another;
                    p_top = p_another;
                    Real rho_e = p_top / (gamma_ - 1.0);
                    vel_top[0] = u_another;
                    vel_top[1] = v_another;
                    E_top = rho_e + 0.5 * rho_top * vel_top.squaredNorm();
                }
                else
                {
                    rho_top = rho0_one;
                    p_top = p_one;
                    Real rho_e = p_top / (gamma_ - 1.0);
                    vel_top[0] = u_one;
                    vel_top[1] = v_one;
                    E_top = rho_e + 0.5 * rho_top * vel_top.squaredNorm();
                }
                CompressibleFluidState state_j(rho_top, vel_top, p_top, E_top);
                CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

                momentum_change_rate -= 2.0 * dW_ijV_j *
				((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
            }
		}
		dmom_dt_[index_i] = momentum_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		mom_[index_i] += dmom_dt_[index_i] * dt;
		vel_[index_i] = mom_[index_i] / rho_[index_i];
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration2ndHalfInFVM<RiemannSolverType>::BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressibleFVM(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_,limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
		Real density_change_rate = 0.0;
		Real energy_change_rate = dE_dt_prior_[index_i];
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			if (inner_neighborhood.boundary_type_[n] == 2)
			{
				CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
				energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
			}

			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				//non-slip wall boundary
				/*Vecd vel_in_wall = -state_i.vel_;
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;*/

				//rigid wall boundary
				Vecd vel_in_wall = (state_i.vel_ - e_ij.dot(state_i.vel_)*(e_ij)) + (-e_ij.dot(state_i.vel_)*(e_ij));
				Real p_in_wall = state_i.p_;
				Real rho_in_wall = state_i.rho_;
				Real E_in_wall=state_i.E_;

				CompressibleFluidState state_j(rho_in_wall, vel_in_wall, p_in_wall, E_in_wall);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
				energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
			}

            if (inner_neighborhood.boundary_type_[n] == 10)
            {
				// parabolic velocity inflow
                //Real U_f = 1.0; //characteristic velocity is set as 1
                //Real h = 1.0;   // the height of the inflow domain size
                //Vecd parabolic_velocity_inlet = Vecd::Zero();
                //parabolic_velocity_inlet[0] = 1.5 * U_f * (1.0 - pos_[index_i][1] * pos_[index_i][1] / pow(0.5 * h, 2));
                //Vecd vel_inlet = parabolic_velocity_inlet;
                //Real p_inlet = state_i.p_;
                //Real rho_inlet = state_i.rho_;

				//given value inlet flow
				Real rho0_another = 8.0;					/**< initial density of another. */
				Real u_another = 8.25*sin(3.14159 / 3.0);	/**< initial velocity of another in X axis. */
				Real v_another = -8.25*cos(3.14159 / 3.0);	/**< initial velocity of another in Y axis. */
				Vecd vel_another= Vecd::Zero();
				vel_another[0] = u_another;
				vel_another[1] = v_another;
				Real p_another = 140.2 / 1.2;					/**< initial pressure of another. */
				Real rho_e_another = p_another / (1.4 - 1.0);
				Real E_inlet_another = rho_e_another + 0.5 * rho0_another * vel_another.squaredNorm();

				Real rho_inlet = rho0_another;
				Real p_inlet = p_another;
				Vecd vel_inlet= Vecd::Zero();
				vel_inlet[0] = u_another;
				vel_inlet[1] = v_another;
				Real E_inlet = E_inlet_another;
                CompressibleFluidState state_j(rho_inlet, vel_inlet, p_inlet,E_inlet);
                CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
				energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
            }
                
            if (inner_neighborhood.boundary_type_[n] == 36)
            {
                //Outlet boundary condition
                Vecd vel_outlet = state_i.vel_;
                Real p_outlet = state_i.p_;
                Real rho_outlet = state_i.rho_;
                Real E_outlet = state_i.E_;
                CompressibleFluidState state_j(rho_outlet, vel_outlet, p_outlet, E_outlet);
				CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

                density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
				energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
            }

			//Top boundary condition
            if (inner_neighborhood.boundary_type_[n] == 4)
            {
                Real rho0_one = 1.4;                         /**< initial density of one fluid. */
                Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
                Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
                Real p_one = 1.0;                            /**< initial pressure of one fluid. */
                Real rho0_another = 8.0;                     /**< initial density of another. */
                Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
                Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
                Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
                Real gamma_ = 1.4;

                Real run_time = GlobalStaticVariables::physical_time_;
                Real x_1 = 1.0 / 6.0 + run_time * 10.0 / sin(3.14159 / 3.0);
                Real p_top = 0.0;
                Real rho_top = 0.0;
                Vecd vel_top = Vecd::Zero();
                Real E_top = 0.0;
                if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - x_1))
                {
                    rho_top = rho0_another;
                    p_top = p_another;
                    Real rho_e = p_top / (gamma_ - 1.0);
                    vel_top[0] = u_another;
                    vel_top[1] = v_another;
                    E_top = rho_e + 0.5 * rho_top * vel_top.squaredNorm();
                }
                else
                {
                    rho_top = rho0_one;
                    p_top = p_one;
                    Real rho_e = p_top / (gamma_ - 1.0);
                    vel_top[0] = u_one;
                    vel_top[1] = v_one;
                    E_top = rho_e + 0.5 * rho_top * vel_top.squaredNorm();
                }
                CompressibleFluidState state_j(rho_top, vel_top, p_top, E_top);
                CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
				energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
            }
		}
		drho_dt_[index_i] = density_change_rate;
		dE_dt_[index_i] = energy_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		E_[index_i] += dE_dt_[index_i] * dt * 0.5;
		rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
	}
	//=================================================================================================//
}
//=================================================================================================//