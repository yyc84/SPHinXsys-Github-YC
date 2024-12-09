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
template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
template <typename... Args>
DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>::
    DiffusionRelaxation(Args &&...args, const StdVec<ContactDiffusionType *> contact_diffusions)
    : DiffusionRelaxation<DataDelegateContact, DiffusionType>(
          std::forward<Args>(args)...), contact_diffusions_(contact_diffusions)
{
    static_assert((... || std::is_same_v < std::decay_t<Args>, DiffusionType),
                  "One of the arguments in args must be of type StdVec<DiffusionType> or DiffusionType");

    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_kernel_gradients_.push_back(ContactKernelGradientType(this->particles_, contact_particles_k));
        contact_Vol_.push_back(contact_particles_k->template registerStateVariable<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
template <typename... Args>
DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>::
    DiffusionRelaxation(Args &&...args, const StdVec<ContactDiffusionType *> contact_diffusions)
    : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>(std::forward<Args>(args)..., contact_diffusions)
{
    contact_gradient_species_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        for (auto &diffusion : this->contact_diffusions_)
        {
            std::string gradient_species_name = diffusion->GradientSpeciesName();
            contact_gradient_species_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
            contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
void DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>::
    getDiffusionChangeRateTwoPhaseHeatExchange(size_t particle_i, size_t particle_j, Vecd &e_ij,
                                    Real surface_area_ij, const StdVec<Real *> &gradient_species_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real thermal_conductivity_i = this->diffusions_[m]->getDiffusionCoeff();
        Real thermal_conductivity_i = this->contact_diffusions_[m]->getDiffusionCoeff();
        Real diff_coeff_ij =
            this->getInterParticleThermalConductivity(thermal_conductivity_i, thermal_conductivity_i);
        Real phi_ij = 2.0 * (this->gradient_species_[m][particle_i] - gradient_species_k[m][particle_j]);
        this->diffusion_dt_[m][particle_i] += diff_coeff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType, class ContactDiffusionType>
void DiffusionRelaxation<TwoPhaseHeatExchange<ContactKernelGradientType>, DiffusionType, ContactDiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &gradient_species_k = this->contact_gradient_species_[k];
        Real *contact_Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij_ = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            getDiffusionChangeRateTwoPhaseHeatExchange(index_i, index_j, e_ij, area_ij, gradient_species_k);
        }
    }
}
//=================================================================================================//
HeatTransferDiffusion::HeatTransferDiffusion(const std::string &diffusion_species_name,
                                       const std::string &gradient_species_name,
                                       Real diff_cf)
    : IsotropicDiffusion(diffusion_species_name, gradient_species_name, diff_cf)
{
    material_type_name_ = "HeatTransferDiffusion";
}
//=================================================================================================//
HeatTransferDiffusion::HeatTransferDiffusion(const std::string &species_name, Real diff_cf)
    : HeatTransferDiffusion(species_name, species_name, diff_cf) {}
//=================================================================================================//
} // namespace SPH
#endif // HEAT_EXCHANGE_TWO_PHASE_HPP