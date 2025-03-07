#include "body_part_by_cell_tracing.h"

#include "base_particles.hpp"
namespace SPH
{
    //=================================================================================================//
    NearShapeSurfaceTracing::
        NearShapeSurfaceTracing(RealBody &real_body, SharedPtr<Shape> shape_ptr, BaseTracingMethod& tracing_cell_method_base)
        : BodyPartByCell(real_body, shape_ptr->getName()), 
        tagging_cell_method_(std::bind(&NearShapeSurfaceTracing::checkNearSurface, this, _1, _2)), tracing_cell_method_base_(tracing_cell_method_base),
      level_set_shape_(level_set_shape_keeper_.createRef<LevelSetShapeLBoundary>(real_body, *shape_ptr.get(), true))
    {
        tagCells(tagging_cell_method_);
    }
    //=================================================================================================//
    bool NearShapeSurfaceTracing::checkNearSurface(Vecd cell_position, Real threshold)
    {
        return level_set_shape_.checkNearSurface(tracing_cell_method_base_.tracingPosition(cell_position, 0.0), threshold);
    }
    //=================================================================================================//
    void NearShapeSurfaceTracing::updateCellList()
    {
        this->body_part_cells_.clear();
        tagCells(tagging_cell_method_);
    }
    //=================================================================================================//
    NearShapeSurfaceStationaryBoundary::NearShapeSurfaceStationaryBoundary(RealBody &real_body, LevelSetShapeLBoundary &level_set_shape)
        : BodyPartByCell(real_body, level_set_shape.getName()), level_set_shape_(level_set_shape)
    {
        TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceStationaryBoundary::checkNearSurface, this, _1, _2);
        tagCells(tagging_cell_method);
    }
    //=================================================================================================//
    NearShapeSurfaceStationaryBoundary::NearShapeSurfaceStationaryBoundary(RealBody &real_body, SharedPtr<Shape> shape_ptr)
        : BodyPartByCell(real_body, shape_ptr->getName()),
          level_set_shape_(level_set_shape_keeper_.createRef<LevelSetShapeLBoundary>(real_body, *shape_ptr.get(), true))
    {
        TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceStationaryBoundary::checkNearSurface, this, _1, _2);
        tagCells(tagging_cell_method);
    }
    //=================================================================================================//
    NearShapeSurfaceStationaryBoundary::NearShapeSurfaceStationaryBoundary(RealBody &real_body)
        : BodyPartByCell(real_body, "NearShapeSurface"),
          level_set_shape_(DynamicCast<LevelSetShapeLBoundary>(this, real_body.getInitialShape()))
    {
        TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceStationaryBoundary::checkNearSurface, this, _1, _2);
        tagCells(tagging_cell_method);
    }
    //=================================================================================================//
    NearShapeSurfaceStationaryBoundary::NearShapeSurfaceStationaryBoundary(RealBody &real_body, const std::string &sub_shape_name)
        : BodyPartByCell(real_body, sub_shape_name),
          level_set_shape_(
              DynamicCast<LevelSetShapeLBoundary>(this, *DynamicCast<ComplexShape>(this, real_body.getInitialShape())
                                                    .getSubShapeByName(sub_shape_name)))
    {
        TaggingCellMethod tagging_cell_method = std::bind(&NearShapeSurfaceStationaryBoundary::checkNearSurface, this, _1, _2);
        tagCells(tagging_cell_method);
    }
    //=================================================================================================//
    bool NearShapeSurfaceStationaryBoundary::checkNearSurface(Vecd cell_position, Real threshold)
    {
        return level_set_shape_.checkNearSurface(cell_position, threshold);
    }
    //=================================================================================================//
} // namespace SPH
