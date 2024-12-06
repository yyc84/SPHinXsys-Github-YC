/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    diffusion_dynamics.h
 * @brief   These are particle dynamics applicable for all type of particles.
 * @author  Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef HEAT_EXCHANGE_TWO_PHASE_H
#define HEAT_EXCHANGE_TWO_PHASE_H

#include "general_diffusion_reaction_dynamics.h"
#include "diffusion_dynamics.h"
#include "diffusion_reaction.h"

namespace SPH
{
/**
 * @class GetDiffusionTimeStepSize
 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
 */


template <typename... InteractionTypes>
class DiffusionRelaxation;


template <class ContactKernelGradientType, class DiffusionType>
class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;
    StdVec<StdVec<Real *>> contact_transfer_;

    void resetContactTransfer(size_t index_i);
    void accumulateDiffusionRate(size_t index_i);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);
    virtual ~DiffusionRelaxation(){};
};

template <typename... ControlTypes>
class Dirichlet; /**< Contact interaction with Dirichlet boundary condition */

template <class ContactKernelGradientType, class DiffusionType>
class DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{

  protected:
    StdVec<StdVec<Real *>> contact_gradient_species_;
    void getDiffusionChangeRateDirichlet(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);
    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

template <typename... ControlTypes>
class TwoPhaseHeatExchange;

template <class ContactKernelGradientType, class DiffusionType>
class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;
    StdVec<StdVec<Real *>> contact_transfer_;

    void resetContactTransfer(size_t index_i);
    void accumulateDiffusionRate(size_t index_i);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);
    virtual ~DiffusionRelaxation(){};
};

template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{

  protected:
    StdVec<DiffusionType *> contact_diffusions_;
    StdVec<StdVec<Real *>> contact_gradient_species_;
    StdVec<StdVec<Real *>> contact_diffusion_dt_;
    StdVec<StdVec<Real *>> contact_thermal_conductivity_;
    void getDiffusionChangeRateTwoPhaseHeatExchange(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);
    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    inline Real getInterParticleThermalConductivity(Real thermal_conductivity_i, Real thermal_conductivity_j)
    {
        return 2 * thermal_conductivity_i * thermal_conductivity_j /( thermal_conductivity_i + thermal_conductivity_j);
    };
};

class HeatTransferDiffusion : public IsotropicDiffusion
{
  protected:
    Real *local_diffusivity_;

  public:
    HeatTransferDiffusion(const std::string &diffusion_species_name,
                            const std::string &gradient_species_name,
                            Real diff_cf = 1.0);
    HeatTransferDiffusion(const std::string &species_name, Real diff_cf = 1.0);
    virtual ~HeatTransferDiffusion(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getReferenceDiffusivity() override { return diff_cf_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return local_diffusivity_[index_i]; };
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
    };
};


} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_H
