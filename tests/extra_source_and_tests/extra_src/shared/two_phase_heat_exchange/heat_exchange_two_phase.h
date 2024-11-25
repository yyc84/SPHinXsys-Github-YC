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

namespace SPH
{
/**
 * @class GetDiffusionTimeStepSize
 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
 */
template <class DiffusionType>
class GetDiffusionTimeStepSize : public BaseDynamics<Real>
{
  public:
    explicit GetDiffusionTimeStepSize(SPHBody &sph_body, DiffusionType &diffusion);
    virtual ~GetDiffusionTimeStepSize(){};

    virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };

  protected:
    Real diff_time_step_;
};

template <typename... InteractionTypes>
class DiffusionRelaxation;

template <class DataDelegationType, class DiffusionType>
class DiffusionRelaxation<DataDelegationType, DiffusionType>
    : public LocalDynamics,
      public DataDelegationType
{
  protected:
    StdVec<DiffusionType *> diffusions_;
    Real *Vol_;
    StdVec<Real *> diffusion_species_;
    StdVec<Real *> gradient_species_;
    StdVec<Real *> diffusion_dt_;

  public:
    template <class BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &body_relation, StdVec<DiffusionType *> diffusions);

    template <class BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &body_relation, DiffusionType *diffusion);

    template <typename BodyRelationType, typename FirstArg>
    explicit DiffusionRelaxation(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : DiffusionRelaxation(parameters.body_relation_, std::get<0>(parameters.others_)){};

    /** So that contact diffusion can be integrated independently without inner interaction. */
    void initialization(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  private:
    void registerSpecies();
};

class KernelGradientInner
{
  public:
    explicit KernelGradientInner(BaseParticles *inner_particles){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return dW_ijV_j * e_ij;
    };
};

class CorrectedKernelGradientInner
{
    Matd *B_;

  public:
    explicit CorrectedKernelGradientInner(BaseParticles *inner_particles)
        : B_(inner_particles->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + B_[index_j]) * e_ij;
    };
};

/**
 * @class DiffusionRelaxationInner
 * @brief Compute the diffusion relaxation process of all species
 */
template <class KernelGradientType, class DiffusionType>
class DiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<DataDelegateInner, DiffusionType>
{
  protected:
    KernelGradientType kernel_gradient_;

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);

    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

class KernelGradientContact
{
  public:
    KernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return dW_ijV_j * e_ij;
    };
};

class CorrectedKernelGradientContact
{
    Matd *B_;
    Matd *contact_B_;

  public:
    CorrectedKernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles)
        : B_(inner_particles->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          contact_B_(contact_particles->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
    };
};

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
class DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{

  protected:
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
};


} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_H
