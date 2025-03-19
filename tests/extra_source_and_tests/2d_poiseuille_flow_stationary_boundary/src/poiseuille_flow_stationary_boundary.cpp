/**
 * @file 	poiseuille_flow.cpp
 * @brief 	2D poiseuille flow example
 * @details This is the one of the basic test cases for validating viscous flow.
 * 			//TODO: this case is too causal now, it should be revised to validate low-Reynolds number flow (Re = 10?).
 * @author 	Chi Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"
#include "level_set_confinement.h"
//#include "io_observation_for_debuging.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 6.0;                   /**< Tank length. */
Real DH = 1.0;                   /**< Tank height. */
Real resolution_ref = DH / 20.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;    /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;
Real mu_f = 1.0e-1;                                    /**< Viscosity. */
Real U_f = 1.0;                                        /**< Reference density of fluid. */
Real gravity_g = 12.0 * mu_f * U_f / rho0_f / DH / DH;            /**< Gravity force of fluid. */                                                /**< Characteristic velocity. */
Real U_max = 1.5 * U_f;                // make sure the maximum anticipated speed
Real c_f = 10.0 * U_max;                                          /**< Reference sound speed. */
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(-2.0 * BW, 0.0));
        inner_wall_shape.push_back(Vecd(-2.0 * BW, DH));
        inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
        inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
        inner_wall_shape.push_back(Vecd(-2.0 * BW, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

class WallUp : public MultiPolygonShape
{
public :
    explicit WallUp(const std::string &shape_name): MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> upper_wall_shape;
        /*upper_wall_shape.push_back(Vecd(-BW, DH + 0.5 * resolution_ref));
        upper_wall_shape.push_back(Vecd(-BW, DH + 1.0 * BW));
        upper_wall_shape.push_back(Vecd(DL + BW, DH + 1.0 * BW));
        upper_wall_shape.push_back(Vecd(DL + BW, DH + 0.5 * resolution_ref));
        upper_wall_shape.push_back(Vecd(-BW, DH + 0.5 * resolution_ref));*/

        upper_wall_shape.push_back(Vecd(-BW, DH ));
        upper_wall_shape.push_back(Vecd(-BW, DH + 1.0 * BW));
        upper_wall_shape.push_back(Vecd(DL + BW, DH + 1.0 * BW));
        upper_wall_shape.push_back(Vecd(DL + BW, DH ));
        upper_wall_shape.push_back(Vecd(-BW, DH ));
         
        multi_polygon_.addAPolygon(upper_wall_shape, ShapeBooleanOps::add);
    }
};

class WallDown : public MultiPolygonShape
{
public :
    explicit WallDown(const std::string &shape_name): MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> down_wall_shape;
        /*down_wall_shape.push_back(Vecd(-BW, - 1.0 * BW));
        down_wall_shape.push_back(Vecd(-BW, 0.0 - 0.5 * resolution_ref));
        down_wall_shape.push_back(Vecd(DL + BW, 0.0 - 0.5 * resolution_ref));
        down_wall_shape.push_back(Vecd(DL + BW, - 1.0 * BW));
        down_wall_shape.push_back(Vecd(-BW, - 1.0 * BW));*/

        down_wall_shape.push_back(Vecd(-BW, -1.0 * BW));
        down_wall_shape.push_back(Vecd(-BW, 0.0));
        down_wall_shape.push_back(Vecd(DL + BW, 0.0 ));
        down_wall_shape.push_back(Vecd(DL + BW, -1.0 * BW));
        down_wall_shape.push_back(Vecd(-BW, -1.0 * BW));

        multi_polygon_.addAPolygon(down_wall_shape, ShapeBooleanOps::add);
    }
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();
    /** topology */
    InnerRelation water_block_inner(water_block);
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    /**
     * @brief 	Methods used for time stepping.
     */
    Gravity gravity(Vecd(gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    
    /**
     * @brief 	Algorithms of fluid dynamics.
     */
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationInner> update_density_by_summation(water_block_inner);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(water_block_inner);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> density_relaxation(water_block_inner);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    /** Impose transport velocity. */
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<TruncatedLinear, AllParticles>> transport_velocity_correction(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<NoLimiter, AllParticles>> transport_velocity_correction(water_block_inner);
    
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> kernel_correction_inner(water_block_inner);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    /** Periodic BCs in x direction. */
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    ParticleSorting particle_sorting(water_block);
    /**
     * @brief Output.
     */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Density");
    ReducedQuantityRecording<TotalKineticEnergy> write_water_kinetic_energy(water_block);

    NearShapeSurfaceStationaryBoundary near_surface_up(water_block, makeShared<InverseShape<WallUp>>("WallUp"));
    near_surface_up.getLevelSetShape().writeLevelSet(sph_system);
    fluid_dynamics::StationaryConfinement confinement_condition_up(near_surface_up);

    NearShapeSurfaceStationaryBoundary near_surface_down(water_block, makeShared<InverseShape<WallDown>>("WallDown"));
    near_surface_down.getLevelSetShape().writeLevelSet(sph_system);
    fluid_dynamics::StationaryConfinement confinement_condition_down(near_surface_down);

    update_density_by_summation.post_processes_.push_back(&confinement_condition_up.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_up.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_up.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_up.surface_bounding_);
    transport_velocity_correction.post_processes_.push_back(&confinement_condition_up.transport_velocity_);
    viscous_force.post_processes_.push_back(&confinement_condition_up.viscous_force_);

    update_density_by_summation.post_processes_.push_back(&confinement_condition_down.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_down.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_down.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_down.surface_bounding_);
    transport_velocity_correction.post_processes_.push_back(&confinement_condition_down.transport_velocity_);
    viscous_force.post_processes_.push_back(&confinement_condition_down.viscous_force_);
    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    constant_gravity.exec();
    /** Output the start states of bodies. */
    body_states_recording.writeToFile(0);
    
    /*ReducedQuantityRecordingForDebuging<Real, ReduceSum<Real>> write_single_variable_real(io_environment, water_block, 0.0, "KernelValueLevelSet");
    GlobalQuantityRecordingForDebuging<Real>wrtie_variable_by_position_real(io_environment, water_block, 0.0, "KernelValueLevelSet");
    ReducedQuantityRecordingForDebuging<Vecd, ReduceSum<Vecd>> write_single_variable_vector(io_environment, water_block, Vecd::Zero(), "KernelGradientLevelSet");
    GlobalQuantityRecordingForDebuging<Vecd>wrtie_variable_by_position_vecd(io_environment, water_block, Vecd::Zero(), "KernelGradientLevelSet");*/
    /**
     * @brief 	Basic parameters.
     */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 100.0;  /**< End time. */
    Real Output_Time = 1.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    /**
     * @brief 	Main loop starts here.
     */
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            /** Acceleration due to viscous force and gravity. */
            time_instance = TickCount::now();

            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            //kernel_correction_inner.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;

            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                write_water_kinetic_energy.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();
            water_block_inner.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        
        /*write_single_variable_real.writeToFile();
        wrtie_variable_by_position_real.writeToFile();
        write_single_variable_vector.writeToFile();
        wrtie_variable_by_position_vecd.writeToFile()*/;
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}
