/**
 * @file 	diffusion_dynamics.hpp
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Y
 */

#ifndef HEAT_EXCHANGE_TWO_PHASE_HPP
#define HEAT_EXCHANGE_TWO_PHASE_HPP

#include "heat_exchange_two_phase.h"

namespace SPH
{


//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(std::forward<Args>(args)...)
{
    contact_gradient_species_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        for (auto &diffusion : this->diffusions_)
        {
            std::string gradient_species_name = diffusion->GradientSpeciesName();
            contact_gradient_species_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
            contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateDirichlet(size_t particle_i, size_t particle_j, Vecd &e_ij,
                                    Real surface_area_ij, const StdVec<Real *> &gradient_species_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real diff_coeff_ij =
            this->diffusions_[m]->getInterParticleDiffusionCoeff(particle_i, particle_i, e_ij);
        Real phi_ij = 2.0 * (this->gradient_species_[m][particle_i] - gradient_species_k[m][particle_j]);
        this->diffusion_dt_[m][particle_i] += diff_coeff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &gradient_species_k = this->contact_gradient_species_[k];
        Real *wall_Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij_ = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            getDiffusionChangeRateDirichlet(index_i, index_j, e_ij, area_ij, gradient_species_k);
        }
    }
}

//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
RungeKuttaStep<DiffusionRelaxationType>::RungeKuttaStep(Args &&...args)
    : DiffusionRelaxationType(std::forward<Args>(args)...)
{
    for (auto &diffusion : this->diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_s_.push_back(
            this->particles_->template registerStateVariable<Real>(diffusion_species_name + "Intermediate"));
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
FirstStageRK2<DiffusionRelaxationType>::FirstStageRK2(Args &&...args)
    : RungeKuttaStep<DiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void FirstStageRK2<DiffusionRelaxationType>::initialization(size_t index_i, Real dt)
{
    DiffusionRelaxationType::initialization(index_i, dt);

    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_species_s_[m][index_i] = this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
SecondStageRK2<DiffusionRelaxationType>::SecondStageRK2(Args &&...args)
    : RungeKuttaStep<DiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void SecondStageRK2<DiffusionRelaxationType>::update(size_t index_i, Real dt)
{
    DiffusionRelaxationType::update(index_i, dt);
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_species_[m][index_i] = 0.5 * this->diffusion_species_s_[m][index_i] +
                                               0.5 * this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename FirstArg, typename... OtherArgs>
DiffusionRelaxationRK2<DiffusionRelaxationType>::
    DiffusionRelaxationRK2(FirstArg &first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      rk2_1st_stage_(first_arg, std::forward<OtherArgs>(other_args)...),
      rk2_2nd_stage_(first_arg, std::forward<OtherArgs>(other_args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void DiffusionRelaxationRK2<DiffusionRelaxationType>::exec(Real dt)
{
    rk2_1st_stage_.exec(dt);
    rk2_2nd_stage_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_HPP