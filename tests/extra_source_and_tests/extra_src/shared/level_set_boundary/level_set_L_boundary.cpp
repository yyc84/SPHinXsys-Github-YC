#include "level_set_L_boundary.h"

#include "adaptation.h"
#include "base_kernel.h"

namespace SPH
{
//=================================================================================================//
MultilevelLevelSetLBoundary::MultilevelLevelSetLBoundary(
    BoundingBox tentative_bounds, Real reference_data_spacing, size_t total_levels,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), kernel_(*sph_adaptation.getKernel()), shape_(shape), total_levels_(total_levels)
{
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(0, reference_data_spacing, global_h_ratio, tentative_bounds);

    for (size_t level = 1; level < total_levels_; ++level) {
        reference_data_spacing *= 0.5;  // Halve the data spacing
        global_h_ratio *= 2;            // Double the ratio
        global_h_ratio_vec_.push_back(global_h_ratio);

        initializeLevel(level, reference_data_spacing, global_h_ratio, tentative_bounds);
    }

    clean_interface_L_boundary = makeUnique<CleanInterfaceLBoundary>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
    correct_topology_L_boundary = makeUnique<CorrectTopologyLBoundary>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
}
//=================================================================================================//
MultilevelLevelSetLBoundary::MultilevelLevelSetLBoundary(
    BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data, Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), kernel_(*sph_adaptation.getKernel()), shape_(shape), total_levels_(1)
{
    Real reference_data_spacing = coarse_data->DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(0, reference_data_spacing, global_h_ratio, tentative_bounds, coarse_data);

    clean_interface_L_boundary = makeUnique<CleanInterfaceLBoundary>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
    correct_topology_L_boundary = makeUnique<CorrectTopologyLBoundary>(*mesh_data_set_.back(), kernel_, global_h_ratio_vec_.back());
}
//=================================================================================================//
void MultilevelLevelSetLBoundary::initializeLevel(size_t level, Real reference_data_spacing, Real global_h_ratio, BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data)
{
    mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));

    RegisterMeshVariableLBoundary register_mesh_variable_L_boundary;
    register_mesh_variable_L_boundary.exec(mesh_data_set_[level]);

    if (coarse_data == nullptr) {
        MeshAllDynamics<InitializeDataInACell> initialize_data_in_a_cell(*mesh_data_set_[level], shape_);
        initialize_data_in_a_cell.exec();
    } else {
        MeshAllDynamics<InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(*mesh_data_set_[level], *coarse_data, shape_);
        initialize_data_in_a_cell_from_coarse.exec();
    }

    FinishDataPackagesLBoundary finish_data_packages_L_boundary(*mesh_data_set_[level], shape_, kernel_, global_h_ratio);
    finish_data_packages_L_boundary.exec();

    registerProbesLBoundary(level);
}
//=================================================================================================//
void MultilevelLevelSetLBoundary::registerProbesLBoundary(size_t level)
{
    probe_signed_distance_set_.push_back(
        probe_signed_distance_vector_keeper_
            .template createPtr<ProbeSignedDistance>(*mesh_data_set_[level]));
    probe_normal_direction_set_.push_back(
        probe_normal_direction_vector_keeper_
            .template createPtr<ProbeNormalDirection>(*mesh_data_set_[level]));
    probe_level_set_gradient_set_.push_back(
        probe_level_set_gradient_vector_keeper_
            .template createPtr<ProbeLevelSetGradient>(*mesh_data_set_[level]));
    probe_kernel_integral_set_.push_back(
        probe_kernel_integral_vector_keeper_
            .template createPtr<ProbeKernelIntegral>(*mesh_data_set_[level]));
    probe_kernel_gradient_integral_set_.push_back(
        probe_kernel_gradient_integral_vector_keeper_
            .template createPtr<ProbeKernelGradientIntegral>(*mesh_data_set_[level]));
    probe_kernel_gradient_multiply_Rij_integral_set_.push_back(
        probe_kernel_multiply_Rij_integral_vector_keeper_
            .template createPtr<ProbeKernelGradientMultiplyRijIntegral>(*mesh_data_set_[level]));
    probe_kernel_gradient_divide_Rij_integral_set_.push_back(
        probe_kernel_divide_Rij_integral_vector_keeper_
            .template createPtr<ProbeKernelGradientDivideRijIntegral>(*mesh_data_set_[level]));
}
//=================================================================================================//
size_t MultilevelLevelSetLBoundary::getCoarseLevel(Real h_ratio)
{
    for (size_t level = total_levels_; level != 0; --level)
        if (h_ratio > global_h_ratio_vec_[level - 1])
            return level - 1; // jump out the loop!

    std::cout << "\n Error: LevelSet level searching out of bound!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    exit(1);
    return 999; // means an error in level searching
};
//=================================================================================================//
void MultilevelLevelSetLBoundary::cleanInterfaceLBoundary(Real small_shift_factor)
{
    clean_interface_L_boundary->exec(small_shift_factor);
}
//=============================================================================================//
void MultilevelLevelSetLBoundary::correctTopologyLBoundary(Real small_shift_factor)
{
    correct_topology_L_boundary->exec(small_shift_factor);
}
//=============================================================================================//
Real MultilevelLevelSetLBoundary::probeSignedDistance(const Vecd &position)
{
    return probe_signed_distance_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
Vecd MultilevelLevelSetLBoundary::probeNormalDirection(const Vecd &position)
{
    return probe_normal_direction_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
Vecd MultilevelLevelSetLBoundary::probeLevelSetGradient(const Vecd &position)
{
    return probe_level_set_gradient_set_[getProbeLevel(position)]->update(position);
}
//=============================================================================================//
size_t MultilevelLevelSetLBoundary::getProbeLevel(const Vecd &position)
{
    for (size_t level = total_levels_; level != 0; --level){
        IsWithinCorePackage is_within_core_package{*mesh_data_set_[level - 1]};
        if(is_within_core_package.update(position))
            return level - 1; // jump out of the loop!
    }
    return 0;
}
//=================================================================================================//
Real MultilevelLevelSetLBoundary::probeKernelIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Real coarse_level_value = probe_kernel_integral_set_[coarse_level]->update(position);
    Real fine_level_value = probe_kernel_integral_set_[coarse_level + 1]->update(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Vecd MultilevelLevelSetLBoundary::probeKernelGradientIntegral(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_gradient_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Vecd coarse_level_value = probe_kernel_gradient_integral_set_[coarse_level]->update(position);
    Vecd fine_level_value = probe_kernel_gradient_integral_set_[coarse_level + 1]->update(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Real MultilevelLevelSetLBoundary::probeKernelGradientMultiplyRij(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_gradient_multiply_Rij_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Real coarse_level_value = probe_kernel_gradient_multiply_Rij_integral_set_[coarse_level]->update(position);
    Real fine_level_value = probe_kernel_gradient_multiply_Rij_integral_set_[coarse_level + 1]->update(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
Real MultilevelLevelSetLBoundary::probeKernelGradientDivideRij(const Vecd &position, Real h_ratio)
{
    if(mesh_data_set_.size() == 1){
        return probe_kernel_gradient_divide_Rij_integral_set_[0]->update(position);
    }
    size_t coarse_level = getCoarseLevel(h_ratio);
    Real alpha = (global_h_ratio_vec_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_vec_[coarse_level + 1] - global_h_ratio_vec_[coarse_level]);
    Real coarse_level_value = probe_kernel_gradient_divide_Rij_integral_set_[coarse_level]->update(position);
    Real fine_level_value = probe_kernel_gradient_divide_Rij_integral_set_[coarse_level + 1]->update(position);

    return alpha * coarse_level_value + (1.0 - alpha) * fine_level_value;
}
//=================================================================================================//
bool MultilevelLevelSetLBoundary::probeIsWithinMeshBound(const Vecd &position)
{
    bool is_bounded = true;
    for (size_t l = 0; l != total_levels_; ++l)
    {
        ProbeIsWithinMeshBound probe_is_within_mesh_bound{*mesh_data_set_[l]};
        if (!probe_is_within_mesh_bound.update(position))
        {
            is_bounded = false;
            break;
        };
    }
    return is_bounded;
}
//=============================================================================================//
} // namespace SPH
