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
 *  HU1527/12-1 and HU1527/12-4												*
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
 * @file 	diffusion_reaction_particles.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef	TWO_PHASE_HEAT_TRANSFER_PARTICLES_H
#define TWO_PHASE_HEAT_TRANSFER_PARTICLES_H

//#include "base_particles.h"
//#include "base_body.h"
//#include "base_material.h"
//#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"
namespace SPH
{

	template <class BaseParticlesType, class DiffusionReactionMaterialType>
	class TwoPhaseDiffusionReactionParticles : public DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>
	{
	public:
		StdLargeVec<Real>  thermal_conductivity_; /**< array of the time derivative of diffusion species */
		StdLargeVec<Real> external_diffusion_dt_sum_;
		StdLargeVec<Real> external_diffusion_dt_;

		TwoPhaseDiffusionReactionParticles(SPHBody& sph_body,
			DiffusionReactionMaterialType* diffusion_reaction_material)
			: DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>(sph_body, diffusion_reaction_material)
		{
			thermal_conductivity_.resize(this->all_species_.size());
			//external_diffusion_dt_sum_.resize(this->all_species_.size());

			external_diffusion_dt_.resize(this->all_species_.size());
			external_diffusion_dt_sum_.resize(this->all_species_.size());

		};
		virtual ~TwoPhaseDiffusionReactionParticles() {};


		virtual void initializeOtherVariables() override
		{
			DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>::initializeOtherVariables();

				//this->registerVariable(external_diffusion_dt_sum_, "HeatFlux");
				this->registerVariable(thermal_conductivity_, "ThermalConductivity");
				//this->template addVariableToWrite<Real>("HeatFlux");
				this->template addVariableToWrite<Real>("ThermalConductivity");
				this->registerVariable(external_diffusion_dt_, "11ExternalChangeRate");
				this->template addVariableToWrite<Real>("11ExternalChangeRate");
				this->registerVariable(external_diffusion_dt_sum_, "HeatFlux");
				this->template addVariableToWrite<Real>("HeatFlux");

		};

		virtual TwoPhaseDiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>* ThisObjectPtr() override { return this; };
	};
}
#endif // TWO_PHASE_HEAT_TRANSFER_PARTICLES_H