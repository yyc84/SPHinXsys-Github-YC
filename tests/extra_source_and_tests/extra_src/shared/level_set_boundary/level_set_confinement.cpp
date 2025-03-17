#include "level_set_confinement.h"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
    //=================================================================================================//
        StationaryConfinementDensity::StationaryConfinementDensity(NearShapeSurfaceStationaryBoundary &near_surface)
        : BaseLocalDynamics<BodyPartByCell>(near_surface),
            rho0_(sph_body_.getBaseMaterial().ReferenceDensity()),
            inv_sigma0_(1.0 / sph_body_.getSPHAdaptation().LatticeNumberDensity()),
            mass_(particles_->getVariableDataByName<Real>("Mass")),
            rho_sum_(particles_->getVariableDataByName<Real>("DensitySummation")),
            pos_(particles_->getVariableDataByName<Vecd>("Position")),
            level_set_shape_(&near_surface.getLevelSetShape()) {}
    //=================================================================================================//
    void StationaryConfinementDensity::update(size_t index_i, Real dt)
    {
        Real inv_Vol_0_i = rho0_ / mass_[index_i];
        rho_sum_[index_i] +=
            level_set_shape_->computeKernelIntegral(pos_[index_i]) * inv_Vol_0_i * rho0_ * inv_sigma0_;
    }
    //=================================================================================================//
    StationaryConfinementIntegration1stHalf::StationaryConfinementIntegration1stHalf(NearShapeSurfaceStationaryBoundary &near_surface)
        : BaseLocalDynamics<BodyPartByCell>(near_surface),
            fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
            rho_(particles_->getVariableDataByName<Real>("Density")),
            p_(particles_->getVariableDataByName<Real>("Pressure")),
            mass_(particles_->getVariableDataByName<Real>("Mass")),
            pos_(particles_->getVariableDataByName<Vecd>("Position")),
            vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
            force_(particles_->getVariableDataByName<Vecd>("Force")),
            level_set_shape_(&near_surface.getLevelSetShape()),
            riemann_solver_(fluid_, fluid_) {}
    //=================================================================================================//
    void StationaryConfinementIntegration1stHalf::update(size_t index_i, Real dt)
    {
        Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
        force_[index_i] -= 2.0 * mass_[index_i] * p_[index_i] * kernel_gradient / rho_[index_i];
    }
    //=================================================================================================//
    StationaryConfinementIntegration2ndHalf::StationaryConfinementIntegration2ndHalf(NearShapeSurfaceStationaryBoundary &near_surface)
        : BaseLocalDynamics<BodyPartByCell>(near_surface),
            fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
            rho_(particles_->getVariableDataByName<Real>("Density")),
            p_(particles_->getVariableDataByName<Real>("Pressure")),
            drho_dt_(particles_->getVariableDataByName<Real>("DensityChangeRate")),
            pos_(particles_->getVariableDataByName<Vecd>("Position")),
            vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
            level_set_shape_(&near_surface.getLevelSetShape()),
            riemann_solver_(fluid_, fluid_) {}
    //=================================================================================================//
    void StationaryConfinementIntegration2ndHalf::update(size_t index_i, Real dt)
    {
        Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
        Vecd vel_j_in_wall = -vel_[index_i];
        drho_dt_[index_i] += rho_[index_i] * (vel_[index_i] - vel_j_in_wall).dot(kernel_gradient);
    }
    //=================================================================================================//
    StationaryConfinementBounding::StationaryConfinementBounding(NearShapeSurfaceStationaryBoundary &near_shape_surface)
        : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          constrained_distance_(0.5 * sph_body_.getSPHAdaptation().MinimumSpacing())
    {
        level_set_shape_ = &near_shape_surface.getLevelSetShape();
    }
    //=================================================================================================//
    void StationaryConfinementBounding::update(size_t index_i, Real dt)
    {
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

        if (phi > -constrained_distance_)
        {
            Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
            pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
        }
    }
	//=================================================================================================//
    template<typename KernelCorrectionType>
    StationaryConfinementTransportVelocity<KernelCorrectionType>
		::StationaryConfinementTransportVelocity(NearShapeSurfaceStationaryBoundary &near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), 
		 fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vecd>("Position")), 
          level_set_shape_(&near_surface.getLevelSetShape()), kernel_correction_(this->particles_),
		zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("ZeroGradientResidue"))
	{}
	//=================================================================================================//
    template <typename KernelCorrectionType>
    void StationaryConfinementTransportVelocity<KernelCorrectionType>::update(size_t index_i, Real dt)
	{
        zero_gradient_residue_[index_i] -= 2.0 * this->kernel_correction_(index_i) *level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
	}
	//=================================================================================================//
    StationaryConfinementTransportVelocitySimple
		::StationaryConfinementTransportVelocitySimple(NearShapeSurfaceStationaryBoundary &near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), 
		 fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vecd>("Position")), 
          level_set_shape_(&near_surface.getLevelSetShape()),
		zero_gradient_residue_(particles_->getVariableDataByName<Vecd>("ZeroGradientResidue"))
	{}
	//=================================================================================================//
    void StationaryConfinementTransportVelocitySimple::update(size_t index_i, Real dt)
	{
        zero_gradient_residue_[index_i] -= 2.0 *level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
	}
	//=================================================================================================//
	//StaticConfinementViscousAcceleration::StaticConfinementViscousAcceleration(NearShapeSurface& near_surface)
	//	: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
	//	pos_(particles_->pos_), mass_(particles_->mass_), force_prior_(particles_->force_prior_), rho_(particles_->rho_),
	//	mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
 //       level_set_shape_(&near_surface.level_set_shape_)
	//	//, force_from_fluid_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForceFromWall")),
	//	//kernel_gradient_rij_(*this->particles_->template registerSharedVariable<Real>("KernelGradientRij"))
	//{
 //       particles_->registerVariable(force_from_fluid_, "ViscousForceFromWall"); 
	//	particles_->registerVariable(kernel_gradient_rij_, "KernelGradientRij");
	//}
	//=================================================================================================//
 //   void StaticConfinementViscousAcceleration::update(size_t index_i, Real dt)
	//{
	//	Vecd force = Vecd::Zero();
	//	Vecd vel_derivative = Vecd::Zero();
	//	Vecd vel_level_set_cell_j = Vecd::Zero();
	//	Real rho_i = rho_[index_i];
	//	/*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
	//	Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
	//	vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
	//	Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
	//	force = 2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
	//	force_prior_[index_i] += force;
 //                   /*below for debuging*/
 //       force_from_fluid_[index_i] = force;
 //       kernel_gradient_rij_[index_i] = kernel_gradient_divide_Rij;
	//	/*for debuging*/
	//	/*Vecd force = Vecd::Zero();
	//	force = 2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative;*/

	//	std::string output_folder = "./output";
	//	std::string filefullpath = output_folder + "/" + "viscous_acceleration_wall_levelset_" + std::to_string(dt) + ".dat";
	//	std::ofstream out_file(filefullpath.c_str(), std::ios::app);
	//	out_file << this->pos_[index_i][0] << " " << this->pos_[index_i][1] << " "<< index_i << " "  << force[0] << " " << force[1]<<" "  << force.norm() << " "
	//	<<kernel_gradient_divide_Rij<< std::endl;

	//}
	//=================================================================================================//
	/*BaseForceFromFluidStaticConfinement::BaseForceFromFluidStaticConfinement(NearShapeSurface& near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_), 
		level_set_shape_(&near_surface.level_set_shape_), force_from_fluid_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForceFromFluid"))
	{
	}*/
	//=================================================================================================//
	//ViscousForceFromFluidStaticConfinement::ViscousForceFromFluidStaticConfinement(NearShapeSurface& near_surface)
	//	:BaseForceFromFluidStaticConfinement(near_surface), pos_(particles_->pos_), rho_(particles_->rho_), mass_(particles_->mass_),
	//	mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_)
	//{
	//	//particles_->registerVariable(force_from_fluid_, "ViscousForceFromFluid");
	//}
	//=================================================================================================//
	/*StaticConfinementExtendIntegration1stHalf::
		StaticConfinementExtendIntegration1stHalf(NearShapeSurface& near_surface, Real penalty_strength)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
		fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), 
		rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
		pos_(particles_->pos_), vel_(particles_->vel_), mass_(particles_->mass_),
		force_(particles_->force_),
		level_set_shape_(&near_surface.level_set_shape_),
		riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}*/
	//=================================================================================================//
	/*void StaticConfinementExtendIntegration1stHalf::update(size_t index_i, Real dt)
	{
		Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
		force_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];
			
		Real penalty_pressure = p_[index_i];
		Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
		Real ratio = distance_to_the_wall  / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
		Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 1.0 * penalty_pressure: 0.0;
			
		force_[index_i] -= 2.0 * penalty_strength_* penalty * kernel_gradient / rho_[index_i];
	}*/
	//=================================================================================================//
	/*StaticConfinementIntegration1stHalfPenaltyVelocity::
		StaticConfinementIntegration1stHalfPenaltyVelocity(NearShapeSurface& near_surface, Real sound_speed, Real penalty_strength)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
		fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), c_0_(sound_speed),
		rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
		pos_(particles_->pos_), vel_(particles_->vel_), mass_(particles_->mass_),
		force_(particles_->force_),
		level_set_shape_(&near_surface.level_set_shape_),
		riemann_solver_(fluid_, fluid_), penalty_strength_(penalty_strength) {}*/
	//=================================================================================================//
	//void StaticConfinementIntegration1stHalfPenaltyVelocity::update(size_t index_i, Real dt)
	//{
	//	Vecd kernel_gradient = level_set_shape_->computeKernelGradientIntegral(pos_[index_i]);
	//	force_[index_i] -= 2.0 * p_[index_i] * kernel_gradient / rho_[index_i];

	//	//Real penalty_pressure = 0.5 * c_0_ * c_0_ * rho_[index_i];
	//	Real penalty_pressure = 0.5 * vel_[index_i].squaredNorm() * rho_[index_i];
	//	Real distance_to_the_wall = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
	//	Real ratio = distance_to_the_wall / (3.0 * sph_body_.sph_adaptation_->ReferenceSpacing());
	//	Real penalty = ratio < 1.0 ? (1.0 - ratio) * (1.0 - ratio) * 0.5 * penalty_pressure : 0.0;

	//	force_[index_i] -= 2.0 * penalty_strength_ * penalty * kernel_gradient / rho_[index_i];
	//}
	//=================================================================================================//
	/*StaticConfinementFreeSurfaceIndication::StaticConfinementFreeSurfaceIndication(NearShapeSurface &near_surface)
	: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
		pos_(particles_->pos_), surface_indicator_(particles_->indicator_),
		level_set_shape_(&near_surface.level_set_shape_), pos_div_(*particles_->getVariableByName<Real>("PositionDivergence"))
	{}*/
	//=================================================================================================//
	//void StaticConfinementFreeSurfaceIndication::interaction(size_t index_i, Real dt )
	//{
	//	//std::string output_folder_1 = "./output";
	//	//std::string filefullpath_1 = output_folder_1 + "/" + "position_divergence_levelset_before" + std::to_string(dt) + ".dat";
	//	//std::ofstream out_file_1(filefullpath_1.c_str(), std::ios::app);
	//	//out_file_1 <<pos_[index_i][0]<<" "<<pos_[index_i][1]<<" " <<index_i<< "  "<<  pos_div_[index_i]<<std::endl;
	//	////out_file << std::fixed << std::setprecision(2)<<pos_[index_i][0]<< "  " <<pos_[index_i][1]<< "  "<< pos_div << "  "<<  pos_div_[index_i]<<std::endl;
	//	//out_file_1 << " \n";

	//	Real pos_div = - level_set_shape_->computeKernelGradientMultiplyRijIntegral(pos_[index_i]);
	//	pos_div_[index_i] += pos_div;
	//	//std::string output_folder = "./output";
	//	//std::string filefullpath = output_folder + "/" + "position_divergence_levelset_after" + std::to_string(dt) + ".dat";
	//	//std::ofstream out_file(filefullpath.c_str(), std::ios::app);
	//	//out_file <<pos_[index_i][0]<<" "<<pos_[index_i][1]<<" " <<index_i<< "  "<<  pos_div_[index_i]<<std::endl;
	//	////out_file << std::fixed << std::setprecision(2)<<pos_[index_i][0]<< "  " <<pos_[index_i][1]<< "  "<< pos_div << "  "<<  pos_div_[index_i]<<std::endl;
	//	//out_file << " \n";
	//}
	//=================================================================================================//
	/*StaticConfinementBounding::StaticConfinementBounding(NearShapeSurface& near_surface)
		: BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
		pos_(particles_->pos_),
		constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
	{
		level_set_shape_ = &near_surface.level_set_shape_;
	}*/
	//=================================================================================================//
	/*void StaticConfinementBounding::update(size_t index_i, Real dt)
	{
		Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

		if (phi > -constrained_distance_)
		{
			Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
			pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
		}
	}*/
	//=================================================================================================//
    StationaryConfinement::StationaryConfinement(NearShapeSurfaceStationaryBoundary &near_surface)
		: density_summation_(near_surface), pressure_relaxation_(near_surface),
          density_relaxation_(near_surface), viscous_force_(near_surface), surface_bounding_(near_surface)
		,transport_velocity_(near_surface)
	{}
	//=================================================================================================//
    StationaryConfinementSimpleMethod::StationaryConfinementSimpleMethod(NearShapeSurfaceStationaryBoundary &near_surface)
		: density_summation_(near_surface), pressure_relaxation_(near_surface),
          density_relaxation_(near_surface), viscous_force_(near_surface), surface_bounding_(near_surface),
          transport_velocity_(near_surface)
	{}
	//=================================================================================================//
   /* StaticConfinementWithPenalty::StaticConfinementWithPenalty(NearShapeSurface &near_surface, Real sound_speed, Real penalty_strength)
		: density_summation_(near_surface), pressure_relaxation_(near_surface),
		density_relaxation_(near_surface), transport_velocity_(near_surface),
		viscous_acceleration_(near_surface), extend_intergration_1st_half_(near_surface, penalty_strength),
		surface_bounding_(near_surface), extend_intergration_1st_half_Velocity(near_surface, sound_speed, penalty_strength)
	{}*/
    //=================================================================================================//
    /*StaticConfinementGeneral::StaticConfinementGeneral(NearShapeSurface &near_surface)
        : density_summation_(near_surface), pressure_relaxation_(near_surface),
            density_relaxation_(near_surface), transport_velocity_(near_surface),
            viscous_acceleration_(near_surface), surface_bounding_(near_surface), free_surface_indication_(near_surface)
    {}*/
		
	}
	//=================================================================================================//
}
//=================================================================================================//