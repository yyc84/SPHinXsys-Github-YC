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

template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactDiffusionType *> contact_diffusions_;
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;

  public:
    template <typename BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &contact_body_relation, DiffusionType *diffusion, StdVec<ContactDiffusionType *> contact_diffusions);

    template <typename BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &body_relation, DiffusionType *diffusion, ContactDiffusionType *contact_diffusion);

    virtual ~DiffusionRelaxation() {};
};

template <typename... ControlTypes>
class HeatExchange;

template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<HeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
{
  protected:
    StdVec<StdVec<Real *>> contact_gradient_species_;

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&...args);

    virtual ~DiffusionRelaxation() {}

    void getDiffusionChangeRateTwoPhaseHeatExchange(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);

    void interaction(size_t index_i, Real dt = 0.0);

    inline Real getInterParticleThermalConductivity(Real thermal_conductivity_i, Real thermal_conductivity_j)
    {
        return 2 * thermal_conductivity_i * thermal_conductivity_j / (thermal_conductivity_i + thermal_conductivity_j);
    };
};

//class HeatTransferDiffusion : public IsotropicDiffusion
//{
//
//  public:
//    HeatTransferDiffusion(const std::string &diffusion_species_name,
//                            const std::string &gradient_species_name,
//                            Real diff_cf = 1.0);
//    HeatTransferDiffusion(const std::string &species_name, Real diff_cf = 1.0);
//    virtual ~HeatTransferDiffusion(){};
//
//    //virtual Real getReferenceDiffusivity() override { return diff_cf_; };
//    //virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return diff_cf_; };
//    /*virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
//    {
//        return diff_cf_;
//    };*/
//    Real getDiffusionCoeff() { return diff_cf_; };
//};

template <typename KernelGradientType, typename ContactKernelGradientType, typename DiffusionType, typename ContactDiffusionType>
class HeatExchangeDiffusionComplex : public LocalDynamics
{
  public:
    using InnerHeatExchange = DiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>;
    using ContactHeatExchange = DiffusionRelaxation<HeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>;

    template <typename BodyRelationInnerType, typename BodyRelationContactType>
    explicit HeatExchangeDiffusionComplex(BodyRelationInnerType &body_relation, BodyRelationContactType &contact_body_relation, DiffusionType *diffusion, ContactDiffusionType *contact_diffusion)
        : HeatExchangeDiffusionComplex(body_relation, contact_body_relation, diffusion, StdVec<ContactDiffusionType *>{contact_diffusion}){};

    template <typename BodyRelationInnerType, typename BodyRelationContactType>
    explicit HeatExchangeDiffusionComplex(BodyRelationInnerType &body_relation, BodyRelationContactType &contact_body_relation, DiffusionType *diffusion, StdVec<ContactDiffusionType *> contact_diffusions)
        : LocalDynamics(body_relation.getSPHBody()), inner_relaxation_(body_relation, diffusion),
          heat_exchange_relaxation_(contact_body_relation, diffusion, contact_diffusions){};

    virtual ~HeatExchangeDiffusionComplex() {};

    void initialization(size_t index_i, Real dt = 0.0)
    {
        inner_relaxation_.initialization(index_i, dt);
    };

    void interaction(size_t index_i, Real dt = 0.0)
    {
        inner_relaxation_.interaction(index_i, dt);
        heat_exchange_relaxation_.interaction(index_i, dt);
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        inner_relaxation_.update(index_i, dt);
    };

  private:
    InnerHeatExchange inner_relaxation_;
    ContactHeatExchange heat_exchange_relaxation_;
};

} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_H
