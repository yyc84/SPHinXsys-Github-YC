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

//template <class InteractionType, class DiffusionType, class ExtraType = void>
//class DiffusionRelaxation;
//
//template <class DataDelegationType, class DiffusionType, class ExtraType = void>
//class DiffusionRelaxation<DataDelegationType, DiffusionType>
//    : public LocalDynamics,
//      public DataDelegationType
//{
//  protected:
//    StdVec<DiffusionType *> diffusions_;
//    Real *Vol_;
//    StdVec<Real *> diffusion_species_;
//    StdVec<Real *> gradient_species_;
//    StdVec<Real *> diffusion_dt_;
//
//  public:
//    template <class BodyRelationType>
//    explicit DiffusionRelaxation(BodyRelationType &body_relation, StdVec<DiffusionType *> diffusions);
//
//    template <class BodyRelationType>
//    explicit DiffusionRelaxation(BodyRelationType &body_relation, DiffusionType *diffusion);
//
//    template <typename BodyRelationType, typename FirstArg>
//    explicit DiffusionRelaxation(ConstructorArgs<BodyRelationType, FirstArg> parameters)
//        : DiffusionRelaxation(parameters.body_relation_, std::get<0>(parameters.others_)){};
//
//    /** So that contact diffusion can be integrated independently without inner interaction. */
//    void initialization(size_t index_i, Real dt = 0.0);
//    void update(size_t index_i, Real dt = 0.0);
//
//  private:
//    void registerSpecies();
//};

//template <typename... InteractionTypes>
//class DiffusionRelaxation;
//
//
//template <class ContactKernelGradientType, class DiffusionType>
//class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
//    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
//{
//  protected:
//    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
//    StdVec<Real *> contact_Vol_;
//    StdVec<StdVec<Real *>> contact_transfer_;
//
//    void resetContactTransfer(size_t index_i);
//    void accumulateDiffusionRate(size_t index_i);
//
//  public:
//    template <typename... Args>
//    explicit DiffusionRelaxation(Args &&... args);
//    virtual ~DiffusionRelaxation(){};
//};

template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactDiffusionType *> contact_diffusions_;
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;

    StdVec<StdVec<Real *>> contact_gradient_species_;
    //StdVec<StdVec<Real *>> contact_transfer_;

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&...args, const StdVec<ContactDiffusionType *> &contact_diffusions)
        : DiffusionRelaxation<DataDelegateContact, DiffusionType>(
          std::forward<Args>(args)...), contact_diffusions_(contact_diffusions)
            {
                for (size_t k = 0; k != this->contact_particles_.size(); ++k)
                {
                    BaseParticles *contact_particles_k = this->contact_particles_[k];
                    contact_kernel_gradients_.push_back(ContactKernelGradientType(this->particles_, contact_particles_k));
                    contact_Vol_.push_back(contact_particles_k->template registerStateVariable<Real>("VolumetricMeasure"));
                }

                contact_gradient_species_.resize(this->contact_particles_.size());
                for (size_t k = 0; k != this->contact_particles_.size(); ++k)
                {
                    BaseParticles *contact_particles_k = this->contact_particles_[k];
                    for (auto &contact_diffusion : this->contact_diffusions_)
                    {
                        std::string gradient_species_name = contact_diffusion->GradientSpeciesName();
                        contact_gradient_species_[k].push_back(
                            contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
                        contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
                    }
                }
            };
          
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&...args, const ContactDiffusionType *contact_diffusion)
        : DiffusionRelaxation(std::forward<Args>(args)..., StdVec<ContactDiffusionType *>{contact_diffusion}){};

    template <typename... Args>
    explicit DiffusionRelaxation(Args &&...args)
        :DiffusionRelaxation<DataDelegateContact, DiffusionType>(std::forward<Args>(args)...){};
    virtual ~DiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0);
    inline Real getInterParticleThermalConductivity(Real thermal_conductivity_i, Real thermal_conductivity_j)
    {
        return 2 * thermal_conductivity_i * thermal_conductivity_j /( thermal_conductivity_i + thermal_conductivity_j);
    };
    void getDiffusionChangeRateTwoPhaseHeatExchange(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);
};


template <typename... ControlTypes>
class TwoPhaseHeatExchange;


template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
{

  protected:
    //StdVec<DiffusionType *> contact_diffusions_;
    StdVec<StdVec<Real *>> contact_gradient_species_;
    //StdVec<StdVec<Real *>> contact_diffusion_dt_;
    //StdVec<StdVec<Real *>> contact_thermal_conductivity_;
    void getDiffusionChangeRateTwoPhaseHeatExchange(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<Real *> &gradient_species_k);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args, const StdVec<ContactDiffusionType *> contact_diffusions)
        : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>(std::forward<Args>(args)...), contact_diffusions_(contact_diffusions)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            BaseParticles *contact_particles_k = this->contact_particles_[k];
            contact_kernel_gradients_.push_back(ContactKernelGradientType(this->particles_, contact_particles_k));
            contact_Vol_.push_back(contact_particles_k->template registerStateVariable<Real>("VolumetricMeasure"));
        }

        contact_gradient_species_.resize(this->contact_particles_.size());
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            BaseParticles *contact_particles_k = this->contact_particles_[k];
            for (auto &contact_diffusion : this->contact_diffusions_)
            {
                std::string gradient_species_name = contact_diffusion->GradientSpeciesName();
                contact_gradient_species_[k].push_back(
                    contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
                contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
            }
        }
    };

    template <typename... Args>
    explicit DiffusionRelaxation(Args &&...args, const ContactDiffusionType * contact_diffusion)
        :DiffusionRelaxation(std::forward<Args>(args)..., StdVec<ContactDiffusionType *>{contact_diffusion}){};
    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    inline Real getInterParticleThermalConductivity(Real thermal_conductivity_i, Real thermal_conductivity_j)
    {
        return 2 * thermal_conductivity_i * thermal_conductivity_j /( thermal_conductivity_i + thermal_conductivity_j);
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

//template <class DiffusionType, class ContactDiffusionType, class KernelGradientType, class ContactKernelGradientType,
//          template <typename... Parameters> typename... ContactInteractionTypes>
//class TwoPhaseHeatExchangeBodyRelaxationComplex
//    : public DiffusionRelaxationRK2<
//          ComplexInteraction<DiffusionRelaxation<
//                                 Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
//                             DiffusionType, ContactDiffusionType>>
//{
//  public:
//    template <typename FirstArg, typename... OtherArgs>
//    explicit TwoPhaseHeatExchangeBodyRelaxationComplex(FirstArg &&first_arg, OtherArgs &&... other_args)
//        : DiffusionRelaxationRK2<
//              ComplexInteraction<DiffusionRelaxation<
//                                     Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
//                                 DiffusionType, ContactDiffusionType>>(first_arg, std::forward<OtherArgs>(other_args)...){};
//    virtual ~TwoPhaseHeatExchangeBodyRelaxationComplex(){};
//};

} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_H
