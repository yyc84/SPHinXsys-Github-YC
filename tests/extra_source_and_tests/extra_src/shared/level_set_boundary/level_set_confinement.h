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
 * @file 	level_set_confinement.h
 * @brief 	Here, we define the Stationary and moving confinement boundary condition classes for fluid dynamics.
 * @details     This boundary condition is based on Level-set filed.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef LEVEL_SET_COFINEMENT_H
#define LEVEL_SET_COFINEMENT_H


#include "shape_confinement.h"
#include "body_part_by_cell_tracing.h"
#include "base_fluid_dynamics.h"
#include "general_constraint.h"
#include "riemann_solver.h"
#include <mutex>
#include "level_set_shape_L_boundary.h"
#include "viscous_dynamics.h"
namespace SPH
{
    namespace fluid_dynamics
    {
    /**
     * @class StationaryConfinementDensity
     * @brief Stationary confinement condition for density summation
     */
    class StationaryConfinementDensity : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        StationaryConfinementDensity(NearShapeSurfaceStationaryBoundary &near_surface);
        virtual ~StationaryConfinementDensity(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real rho0_, inv_sigma0_;
        Real *mass_, *rho_sum_;
        Vecd *pos_;
        LevelSetShapeLBoundary *level_set_shape_;

        /*below for debuging*/
        Real *kernel_weight_ij_, *kernel_weight_wall_ij_;
    };

    /**
     * @class StationaryConfinementIntegration1stHalf
     * @brief Stationary confinement condition for pressure relaxation
     */
    class StationaryConfinementIntegration1stHalf : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        StationaryConfinementIntegration1stHalf(NearShapeSurfaceStationaryBoundary &near_surface);
        virtual ~StationaryConfinementIntegration1stHalf(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Fluid &fluid_;
        Real *rho_, *p_, *mass_;
        Vecd *pos_, *vel_, *force_;
        LevelSetShapeLBoundary *level_set_shape_;
        AcousticRiemannSolver riemann_solver_;

        /*below for debuging*/ 
        Vecd *kernel_gradient_ij_, *kernel_gradient_wall_ij_;
    };

    /**
     * @class StationaryConfinementIntegration2ndHalf
     * @brief Stationary confinement condition for density relaxation
     */
    class StationaryConfinementIntegration2ndHalf : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        StationaryConfinementIntegration2ndHalf(NearShapeSurfaceStationaryBoundary &near_surface);
        virtual ~StationaryConfinementIntegration2ndHalf(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Fluid &fluid_;
        Real *rho_, *p_, *drho_dt_;
        Vecd *pos_, *vel_;
        LevelSetShapeLBoundary *level_set_shape_;
        AcousticRiemannSolver riemann_solver_;
    };

     /**
     * @class StationaryConfinementBounding
     * @brief Stationary confinement condition for pressure relaxation
     * map constrained particles to level-set boundary and
     * r = r + phi * norm (vector distance to face)
     */
    class StationaryConfinementBounding : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        StationaryConfinementBounding(NearShapeSurfaceStationaryBoundary &near_surface);
        virtual ~StationaryConfinementBounding(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *pos_;
        LevelSetShapeLBoundary *level_set_shape_;
        Real constrained_distance_;
    };

        /**
         * @class StationaryConfinementTransportVelocity
         * @brief Stationary confinement condition for transport velocity
         */
    template <typename KernelCorrectionType>
    class StationaryConfinementTransportVelocity : public BaseLocalDynamics<BodyPartByCell>
        {
        public:
            StationaryConfinementTransportVelocity(NearShapeSurfaceStationaryBoundary &near_surface);
          virtual ~StationaryConfinementTransportVelocity(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Fluid &fluid_;
            Vecd *pos_;
            //int *surface_indicator_;
            //const Real coefficient_;
            //Real smoothing_length_sqr_;
            
            LevelSetShapeLBoundary *level_set_shape_;
            Vecd *zero_gradient_residue_;
            KernelCorrectionType kernel_correction_;
        };

    using StationaryConfinementTransportVelocityLinearCorrection = StationaryConfinementTransportVelocity<LinearGradientCorrection>; 

        /**
         * @class StationaryConfinementViscousAcceleration
         * @brief Stationary confinement condition for viscous acceleration
         */
        //class StationaryConfinementViscousAcceleration : public BaseLocalDynamics<BodyPartByCell>
        //{
        //public:
        //    StationaryConfinementViscousAcceleration(NearShapeSurfaceStationaryBoundary &near_surface);
        //  virtual ~StationaryConfinementViscousAcceleration(){};
        //    void update(size_t index_i, Real dt = 0.0);
        //    //StdLargeVec<Vecd> &getForceFromFluid() { return force_from_fluid_; };
        //protected:
        //    Vecd *pos_;
        //    Real *rho_, *mass_;
        //    Vecd *vel_, *force_prior_;
        //    Real mu_;
        //    LevelSetShape* level_set_shape_;
        //    Vecd *force_from_fluid_;
        //    Real *kernel_gradient_rij_;

        //};

     class StationaryConfinementTransportVelocitySimple : public BaseLocalDynamics<BodyPartByCell>
    {
      public:
        StationaryConfinementTransportVelocitySimple(NearShapeSurfaceStationaryBoundary &near_surface);
        virtual ~StationaryConfinementTransportVelocitySimple(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Fluid &fluid_;
        Vecd *pos_;
        // int *surface_indicator_;
        // const Real coefficient_;
        // Real smoothing_length_sqr_;

        LevelSetShapeLBoundary *level_set_shape_;
        Vecd *zero_gradient_residue_;
        //KernelCorrectionType kernel_correction_;
    };

        /**
         * @class StationaryConfinementViscousAcceleration
         * @brief Stationary confinement condition for viscous acceleration
         */
        template <typename ViscosityType>
        class StationaryConfinementViscousForce : public BaseLocalDynamics<BodyPartByCell>
        {
        public:
            StationaryConfinementViscousForce(NearShapeSurfaceStationaryBoundary &near_surface)
                : BaseLocalDynamics<BodyPartByCell>(near_surface),
                pos_(particles_->getVariableDataByName<Vecd>("Position")),
                mass_(particles_->getVariableDataByName<Real>("Mass")), 
                viscous_force_(particles_->getVariableDataByName<Vecd>("ViscousForce")),
                rho_(particles_->getVariableDataByName<Real>("Density")),
                mu_(particles_),
                vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
                fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), viscous_force_on_solid_(particles_->template registerStateVariable<Vecd>("ViscousForceOnWall")),
                viscous_force_from_fluid_(particles_->template registerStateVariable<Vecd>("ViscousForceFromFluid")),
                level_set_shape_(&near_surface.getLevelSetShape())
		{}
            virtual ~StationaryConfinementViscousForce(){};
            void update(size_t index_i, Real dt = 0.0)
            {
                Vecd force = Vecd::Zero();
			    Vecd vel_derivative = Vecd::Zero();
			    Vecd vel_level_set_cell_j = Vecd::Zero();
			    Real rho_i = rho_[index_i];
			    /*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			    //Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			    vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			    Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
                            force = 2.0 * mu_(index_i, index_i) * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
                viscous_force_[index_i] += force;
                viscous_force_on_solid_[index_i] = force;
                viscous_force_from_fluid_[index_i] = -force;
            }
            //Vecd& getForceFromFluid() { return viscous_force_; };
        protected:
            Fluid &fluid_;
            Vecd *pos_, *vel_, *viscous_force_;
            Real *rho_, *mass_;
            
            ViscosityType mu_;
            LevelSetShapeLBoundary* level_set_shape_;
            //Vecd *force_from_fluid_;
            Vecd *viscous_force_on_solid_;
            Vecd *viscous_force_from_fluid_;
        };

        using StationaryConfinementFixedViscousForce = StationaryConfinementViscousForce<FixedViscosity>;

        class StationaryConfinementViscousForceSimple : public BaseLocalDynamics<BodyPartByCell>
        {
          public:
            StationaryConfinementViscousForceSimple(NearShapeSurfaceStationaryBoundary &near_surface )
                : BaseLocalDynamics<BodyPartByCell>(near_surface),
                  pos_(particles_->getVariableDataByName<Vecd>("Position")),
                  mass_(particles_->getVariableDataByName<Real>("Mass")),
                  viscous_force_(particles_->getVariableDataByName<Vecd>("ViscousForce")),
                  rho_(particles_->getVariableDataByName<Real>("Density")),
                  mu_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial()).ReferenceViscosity()),
                  vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
                  fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())), viscous_force_on_solid_(particles_->template registerStateVariable<Vecd>("ViscousForceOnWall")),
                  viscous_force_from_fluid_(particles_->template registerStateVariable<Vecd>("ViscousForceFromFluid")),
                  level_set_shape_(&near_surface.getLevelSetShape())
            {
            }
            virtual ~StationaryConfinementViscousForceSimple(){};
            void update(size_t index_i, Real dt = 0.0)
            {
                Vecd force = Vecd::Zero();
                Vecd vel_derivative = Vecd::Zero();
                Vecd vel_level_set_cell_j = Vecd::Zero();
                Real rho_i = rho_[index_i];
                /*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
                // Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
                vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
                Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
                force = 2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative / rho_i;
                viscous_force_[index_i] += force;
                viscous_force_on_solid_[index_i] = force;
                viscous_force_from_fluid_[index_i] = -force;
            }
            // Vecd& getForceFromFluid() { return viscous_force_; };
          protected:
            Fluid &fluid_;
            Vecd *pos_, *vel_, *viscous_force_;
            Real *rho_, *mass_;

            Real mu_;
            LevelSetShapeLBoundary *level_set_shape_;
            Vecd *viscous_force_on_solid_;
            Vecd *viscous_force_from_fluid_;
        };

        //class BaseForceFromFluidStationaryConfinement : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        //{
        //public: 
        //    explicit BaseForceFromFluidStationaryConfinement(NearShapeSurface &near_surface);
        //  virtual ~BaseForceFromFluidStationaryConfinement(){};
        //    StdLargeVec<Vecd> &getForceFromFluid() { return force_from_fluid_; };

        //protected:
        //    StdLargeVec<Vecd> force_from_fluid_;
        //    LevelSetShape* level_set_shape_;
        //    //StdLargeVec<Real> &Vol_;
        //};

       // class ViscousForceFromFluidStationaryConfinement : public BaseForceFromFluidStationaryConfinement
       // {
       // public:
       //     explicit ViscousForceFromFluidStationaryConfinement(NearShapeSurface &near_surface);
       //   virtual ~ViscousForceFromFluidStationaryConfinement(){};
       //     inline void interaction(size_t index_i, Real dt = 0.0)
       //     {
			    //Vecd acceleration = Vecd::Zero();
			    //Vecd vel_derivative = Vecd::Zero();
			    //Vecd vel_level_set_cell_j = Vecd::Zero();
			    //Real rho_i = rho_[index_i];
			    ///*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			    //vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			    //Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
			    //force_from_fluid_[index_i]= -2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
       //         //force_from_fluid_[index_i]= -2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative ;
       //     }
       // protected:
       //     StdLargeVec<Vecd>& pos_;
       //     StdLargeVec<Real>& rho_, &mass_;
       //     StdLargeVec<Vecd>& vel_;
       //     Real mu_;
       // };
        /**
        * @class StationaryConfinementIntegration1stHalf
        * @brief Stationary confinement condition for pressure relaxation
        */
        /*class StationaryConfinementExtendIntegration1stHalf : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StationaryConfinementExtendIntegration1stHalf(NearShapeSurface &near_surface, Real penalty_strength = 2.0);
          virtual ~StationaryConfinementExtendIntegration1stHalf(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real penalty_strength_;
            Fluid& fluid_;
            StdLargeVec<Real>&rho_, &p_, &mass_;
            StdLargeVec<Vecd>&pos_, &vel_, &force_;
            LevelSetShape* level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };*/

        /**
        * @class StaticConfinementIntegration1stHalf
        * @brief static confinement condition for pressure relaxation
        */
        /*class StationaryConfinementIntegration1stHalfPenaltyVelocity : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StationaryConfinementIntegration1stHalfPenaltyVelocity(NearShapeSurface &near_surface, Real sound_speed, Real penalty_strength = 2.0);
          virtual ~StationaryConfinementIntegration1stHalfPenaltyVelocity(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real penalty_strength_, c_0_;
            Fluid& fluid_;
            StdLargeVec<Real>& rho_, & p_, &mass_;
            StdLargeVec<Vecd>& pos_, & vel_, &force_;
            LevelSetShape* level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };*/

         /**
        * @class StationaryConfinementFreeSurfaceIndication
        * @brief Stationary confinement condition for free surface particle indicate
        */
        /*class StationaryConfinementFreeSurfaceIndication : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StationaryConfinementFreeSurfaceIndication(NearShapeSurface &near_surface);
          virtual ~StationaryConfinementFreeSurfaceIndication(){};
            void interaction(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<Real>& pos_div_;
            StdLargeVec<int> &surface_indicator_;
            LevelSetShape* level_set_shape_;
        };*/

       

        /**
        * @class StaticConfinement
        * @brief Static confined boundary condition for complex structures with bounding.
        */
        class StationaryConfinement
         {
         public:

           SimpleDynamics<StationaryConfinementDensity> density_summation_;
           SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
           SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
           SimpleDynamics<StationaryConfinementBounding> surface_bounding_;
           SimpleDynamics<StationaryConfinementFixedViscousForce> viscous_force_;
           SimpleDynamics<StationaryConfinementTransportVelocityLinearCorrection> transport_velocity_;

         StationaryConfinement(NearShapeSurfaceStationaryBoundary &near_surface);
           virtual ~StationaryConfinement(){};
         };

        class StationaryConfinementWithOutTransportVelocit
         {
           public:
             SimpleDynamics<StationaryConfinementDensity> density_summation_;
             SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
             SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
             SimpleDynamics<StationaryConfinementBounding> surface_bounding_;
             SimpleDynamics<StationaryConfinementFixedViscousForce> viscous_force_;

             StationaryConfinementWithOutTransportVelocit(NearShapeSurfaceStationaryBoundary &near_surface):
                 density_summation_(near_surface), pressure_relaxation_(near_surface),
                 density_relaxation_(near_surface), viscous_force_(near_surface), surface_bounding_(near_surface){}
             virtual ~StationaryConfinementWithOutTransportVelocit(){};
         };

        class StationaryConfinementWithOutViscousForce
         {
           public:
             SimpleDynamics<StationaryConfinementDensity> density_summation_;
             SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
             SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
             SimpleDynamics<StationaryConfinementBounding> surface_bounding_;
             SimpleDynamics<StationaryConfinementTransportVelocitySimple> transport_velocity_;

             StationaryConfinementWithOutViscousForce(NearShapeSurfaceStationaryBoundary &near_surface) : 
                 density_summation_(near_surface), pressure_relaxation_(near_surface),
                density_relaxation_(near_surface), transport_velocity_(near_surface), surface_bounding_(near_surface) {}
             virtual ~StationaryConfinementWithOutViscousForce(){};
         };

        class StationaryConfinementSimple
         {
           public:
             SimpleDynamics<StationaryConfinementDensity> density_summation_;
             SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
             SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
             SimpleDynamics<StationaryConfinementBounding> surface_bounding_;

             StationaryConfinementSimple(NearShapeSurfaceStationaryBoundary &near_surface) : 
                 density_summation_(near_surface), pressure_relaxation_(near_surface),
                  density_relaxation_(near_surface),  surface_bounding_(near_surface) {}
             virtual ~StationaryConfinementSimple(){};
         };

        class StationaryConfinementSimpleMethod
         {
           public:
             SimpleDynamics<StationaryConfinementDensity> density_summation_;
             SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
             SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
             SimpleDynamics<StationaryConfinementBounding> surface_bounding_;
             SimpleDynamics<StationaryConfinementViscousForceSimple> viscous_force_;
             SimpleDynamics<StationaryConfinementTransportVelocitySimple> transport_velocity_;

             StationaryConfinementSimpleMethod(NearShapeSurfaceStationaryBoundary &near_surface);
             virtual ~StationaryConfinementSimpleMethod(){};
         };

        /**
        * @class StationaryConfinement
        * @brief Stationary confined boundary condition for complex structures with penalty force for light phase.
        */
        /* class StationaryConfinementWithPenalty
        {
        public:

            SimpleDynamics<StationaryConfinementDensity> density_summation_;
          SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
            SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
          SimpleDynamics<StationaryConfinementTransportVelocity> transport_velocity_;
            SimpleDynamics<StationaryConfinementViscousAcceleration, SequencedPolicy> viscous_acceleration_;
          SimpleDynamics<StationaryConfinementExtendIntegration1stHalf> extend_intergration_1st_half_;
            SimpleDynamics<StationaryConfinementIntegration1stHalfPenaltyVelocity> extend_intergration_1st_half_Velocity;
          SimpleDynamics<StationaryConfinementBounding> surface_bounding_;

            StationaryConfinementWithPenalty(NearShapeSurface &near_surface, Real sound_speed, Real penalty_strength);
          virtual ~StationaryConfinementWithPenalty(){};
        };

        class StationaryConfinementGeneral
        {
        public:
            SimpleDynamics<StationaryConfinementDensity, SequencedPolicy> density_summation_;
          SimpleDynamics<StationaryConfinementIntegration1stHalf> pressure_relaxation_;
            SimpleDynamics<StationaryConfinementIntegration2ndHalf> density_relaxation_;
          SimpleDynamics<StationaryConfinementTransportVelocity, SequencedPolicy> transport_velocity_;
            SimpleDynamics<StationaryConfinementViscousAcceleration, SequencedPolicy> viscous_acceleration_;
          InteractionDynamics<StationaryConfinementFreeSurfaceIndication> free_surface_indication_;
            SimpleDynamics<StationaryConfinementBounding> surface_bounding_;

            StationaryConfinementGeneral(NearShapeSurface &near_surface);
            virtual ~StationaryConfinementGeneral(){};
        };*/

    }
}
#endif  LEVEL_SET_COFINEMENT_H
