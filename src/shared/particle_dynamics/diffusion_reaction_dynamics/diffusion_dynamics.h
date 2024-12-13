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

#ifndef DIFFUSION_DYNAMICS_H
#define DIFFUSION_DYNAMICS_H

#include "general_diffusion_reaction_dynamics.h"

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

template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
    : public DiffusionRelaxation<DataDelegateContact, DiffusionType>
{
  protected:
    StdVec<ContactDiffusionType *> contact_diffusions_;
    //ContactDiffusionType *contact_diffusion_;
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;
    StdVec<Real *> contact_Vol_;

  public:
    //template <typename... Args>
    template <typename BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &contact_body_relation, DiffusionType *diffusion, StdVec<ContactDiffusionType *> contact_diffusions);
   
    template <typename BodyRelationType>
    explicit DiffusionRelaxation(BodyRelationType &body_relation, DiffusionType *diffusion, ContactDiffusionType *contact_diffusion);

    virtual ~DiffusionRelaxation(){};
};

//template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
//class DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
//    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
//{
//  protected:
//    StdVec<ContactDiffusionType *> contact_diffusions_;
//
//  public:
//    template <typename... Args>
//    explicit DiffusionRelaxation(Args &&...args, DiffusionType &diffusion, const StdVec<ContactDiffusionType *> &contact_diffusions)
//        : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(std::forward<Args>(args)..., diffusion),
//          contact_diffusions_(contact_diffusions) {}
//
//    template <typename... Args>
//    explicit DiffusionRelaxation(Contact<ContactKernelGradientType> &contact, DiffusionType &diffusion, ContactDiffusionType *contact_diffusion)
//        : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(std::forward<Args>(args)..., diffusion),
//          contact_diffusions_({contact_diffusion}) {}
//
//    virtual ~DiffusionRelaxation() {}
//};


//template <typename... ControlTypes>
//class HeatExchange;
//
//template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
//class DiffusionRelaxation<HeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
//    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>
//{
//  protected:
//    StdVec<StdVec<Real *>> contact_gradient_species_;
//
//  public:
//    template <typename BodyRelationType>
//    explicit DiffusionRelaxation(Args &&...args, DiffusionType *diffusion, const StdVec<ContactDiffusionType *> &contact_diffusions)
//        : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>(std::forward<Args>(args)..., diffusion, contact_diffusions)
//    {
//        contact_gradient_species_.resize(this->contact_particles_.size());
//        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
//        {
//            BaseParticles *contact_particles_k = this->contact_particles_[k];
//            for (auto &contact_diffusion : this->contact_diffusions_)
//            {
//                std::string gradient_species_name = contact_diffusion->GradientSpeciesName();
//                contact_gradient_species_[k].push_back(
//                    contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
//                contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
//            };
//        };
//    };
//
//    virtual ~DiffusionRelaxation() {}
//
//    void getDiffusionChangeRateTwoPhaseHeatExchange(
//        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
//        const StdVec<Real *> &gradient_species_k)
//    {
//        for (size_t m = 0; m < this->diffusions_.size(); ++m)
//        {
//            Real thermal_conductivity_i = this->diffusions_[m]->getReferenceDiffusivity();
//            Real thermal_conductivity_j = this->contact_diffusions_[m]->getReferenceDiffusivity();
//            Real diff_coeff_ij =
//                this->getInterParticleThermalConductivity(thermal_conductivity_i, thermal_conductivity_j);
//            Real phi_ij = 2.0 * (this->gradient_species_[m][particle_i] - gradient_species_k[m][particle_j]);
//            this->diffusion_dt_[m][particle_i] += diff_coeff_ij * phi_ij * surface_area_ij;
//        };
//    };
//
//    void interaction(size_t index_i, Real dt = 0.0)
//    {
//        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
//        {
//            StdVec<Real *> &gradient_species_k = this->contact_gradient_species_[k];
//            Real *contact_Vol_k = this->contact_Vol_[k];
//            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
//            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
//            {
//                size_t index_j = contact_neighborhood.j_[n];
//                Real r_ij_ = contact_neighborhood.r_ij_[n];
//                Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_Vol_k[index_j];
//                Vecd &e_ij = contact_neighborhood.e_ij_[n];
//
//                const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
//                Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
//                getDiffusionChangeRateTwoPhaseHeatExchange(index_i, index_j, e_ij, area_ij, gradient_species_k);
//            };
//        };
//    };
//
//    inline Real getInterParticleThermalConductivity(Real thermal_conductivity_i, Real thermal_conductivity_j)
//    {
//        return 2 * thermal_conductivity_i * thermal_conductivity_j / (thermal_conductivity_i + thermal_conductivity_j);
//    };
//};

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
class Neumann; /**< Contact interaction with Neumann boundary condition */

template <class ContactKernelGradientType, class DiffusionType>
class DiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{
    Vecd *n_;
    StdVec<StdVec<Real *>> contact_diffusive_flux_;
    StdVec<Vecd *> contact_n_;

  protected:
    void getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j,
                                       Real surface_area_ij_Neumann,
                                       const StdVec<Real *> &diffusive_flux_k);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);
    virtual ~DiffusionRelaxation(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <typename... ControlTypes>
class Robin; /**< Contact interaction with Robin boundary condition */

template <class ContactKernelGradientType, class DiffusionType>
class DiffusionRelaxation<Robin<ContactKernelGradientType>, DiffusionType>
    : public DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>
{
    Vecd *n_;
    StdVec<StdVec<Real *>> contact_convection_;
    StdVec<StdVec<Real *>> contact_species_infinity_;
    StdVec<Vecd *> contact_n_;

  protected:
    void getTransferRateRobin(
        size_t particle_i, size_t particle_j, Real surface_area_ij_Robin,
        StdVec<Real *> &transfer_k,
        StdVec<Real *> &convection_k,
        StdVec<Real *> &species_infinity_k);

  public:
    template <typename... Args>
    explicit DiffusionRelaxation(Args &&... args);

    virtual ~DiffusionRelaxation(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class RungeKuttaStep
 * @brief A general step for runge-kutta integration scheme.
 * @details Am intermediate state for species is introduced here
 * to achieve multi-step integration.
 */
template <class DiffusionRelaxationType>
class RungeKuttaStep : public DiffusionRelaxationType
{
  protected:
    StdVec<Real *> diffusion_species_s_;

  public:
    template <typename... Args>
    RungeKuttaStep(Args &&... args);

    virtual ~RungeKuttaStep(){};
};

/**
 * @class FirstStageRK2
 * @brief The first stage of a 2nd-order runge-kutta integration scheme.
 * A intermediate state for species is introduced here to achieve multi-step integration.
 */
template <class DiffusionRelaxationType>
class FirstStageRK2 : public RungeKuttaStep<DiffusionRelaxationType>
{
  public:
    template <typename... Args>
    FirstStageRK2(Args &&... args);

    virtual ~FirstStageRK2(){};
    void initialization(size_t index_i, Real dt = 0.0);
};

/**
 * @class SecondStageRK2
 * @brief The second stage of the 2nd-order Runge-Kutta scheme.
 */
template <class DiffusionRelaxationType>
class SecondStageRK2 : public RungeKuttaStep<DiffusionRelaxationType>
{
  public:
    template <typename... Args>
    SecondStageRK2(Args &&... args);

    virtual ~SecondStageRK2(){};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class DiffusionRelaxationRK2
 * @brief The 2nd-order runge-kutta integration scheme.
 */
template <class DiffusionRelaxationType>
class DiffusionRelaxationRK2 : public BaseDynamics<void>
{
  protected:
    Dynamics1Level<FirstStageRK2<DiffusionRelaxationType>, SequencedPolicy> rk2_1st_stage_;
    Dynamics1Level<SecondStageRK2<DiffusionRelaxationType>, SequencedPolicy> rk2_2nd_stage_;

  public:
    template <typename FirstArg, typename... OtherArgs>
    explicit DiffusionRelaxationRK2(FirstArg &first_arg, OtherArgs &&... other_args);

    virtual ~DiffusionRelaxationRK2(){};

    virtual void exec(Real dt = 0.0) override;
};

template <class DiffusionType, class KernelGradientType, class ContactKernelGradientType,
          template <typename... Parameters> typename... ContactInteractionTypes>
class DiffusionBodyRelaxationComplex
    : public DiffusionRelaxationRK2<
          ComplexInteraction<DiffusionRelaxation<
                                 Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
                             DiffusionType>>
{
  public:
    template <typename FirstArg, typename... OtherArgs>
    explicit DiffusionBodyRelaxationComplex(FirstArg &&first_arg, OtherArgs &&... other_args)
        : DiffusionRelaxationRK2<
              ComplexInteraction<DiffusionRelaxation<
                                     Inner<KernelGradientType>, ContactInteractionTypes<ContactKernelGradientType>...>,
                                 DiffusionType>>(first_arg, std::forward<OtherArgs>(other_args)...){};
    virtual ~DiffusionBodyRelaxationComplex(){};
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_H
