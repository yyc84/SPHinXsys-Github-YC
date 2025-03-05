#include "level_set_shape_L_boundary.h"

#include "base_body.h"
#include "io_all.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
LevelSetShapeLBoundary::
    LevelSetShapeLBoundary(Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio)
    : Shape(shape.getName()), sph_adaptation_(sph_adaptation),
      level_set_L_boundary_(*level_set_L_boundary_keeper_.movePtr(sph_adaptation->createLevelSetLBoundary(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
LevelSetShapeLBoundary::LevelSetShapeLBoundary(SPHBody &sph_body, Shape &shape, Real refinement_ratio)
    : Shape(shape.getName()),
      level_set_L_boundary_(*level_set_L_boundary_keeper_.movePtr(
          sph_body.getSPHAdaptation().createLevelSetLBoundary(shape, refinement_ratio)))
{
    bounding_box_ = shape.getBounds();
    is_bounds_found_ = true;
}
//=================================================================================================//
void LevelSetShapeLBoundary::writeLevelSet(SPHSystem &sph_system)
{
    MeshRecordingToPlt write_level_set_to_plt(sph_system, level_set_L_boundary_);
    write_level_set_to_plt.writeToFile(0);
}
//=================================================================================================//
LevelSetShapeLBoundary *LevelSetShapeLBoundary::cleanLevelSet(Real small_shift_factor)
{
    level_set_L_boundary_.cleanInterfaceLBoundary(small_shift_factor);
    return this;
}
//=================================================================================================//
LevelSetShapeLBoundary *LevelSetShapeLBoundary::correctLevelSetSign(Real small_shift_factor)
{
    level_set_L_boundary_.correctTopologyLBoundary(small_shift_factor);
    return this;
}
//=================================================================================================//
bool LevelSetShapeLBoundary::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return level_set_L_boundary_.probeSignedDistance(probe_point) < 0.0 ? true : false;
}
//=================================================================================================//
Vecd LevelSetShapeLBoundary::findClosestPoint(const Vecd &probe_point)
{
    Real phi = level_set_L_boundary_.probeSignedDistance(probe_point);
    Vecd normal = level_set_L_boundary_.probeNormalDirection(probe_point);
    return probe_point - phi * normal;
}
//=================================================================================================//
BoundingBox LevelSetShapeLBoundary::findBounds()
{
    if (!is_bounds_found_)
    {
        std::cout << "\n FAILURE: LevelSetShape bounds should be defined at construction!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return bounding_box_;
}
//=================================================================================================//
Vecd LevelSetShapeLBoundary::findLevelSetGradient(const Vecd &probe_point)
{
    return level_set_L_boundary_.probeLevelSetGradient(probe_point);
}
//=================================================================================================//
Real LevelSetShapeLBoundary::computeKernelIntegral(const Vecd &probe_point, Real h_ratio)
{
    return level_set_L_boundary_.probeKernelIntegral(probe_point, h_ratio);
}
//=================================================================================================//
Vecd LevelSetShapeLBoundary::computeKernelGradientIntegral(const Vecd &probe_point, Real h_ratio)
{
    return level_set_L_boundary_.probeKernelGradientIntegral(probe_point, h_ratio);
}
//=================================================================================================//
} // namespace SPH
