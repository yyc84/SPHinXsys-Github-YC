/**
 * @file 	2d_flow_around_cylinder_static_confinement.cpp
 * @brief 	This is the benchmark test for the wall modeling of viscous flow.
 * @details We consider a flow passing by a cylinder in 2D.
 * @author 	Yongchuan Yu
 */
#include "2d_flow_around_cylinder_static_confinement.h"
#include "sphinxsys.h"
#include "level_set_confinement.h"
//#include "io_observation_for_debuging.h"
using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    ParameterizationIO *parameterization_io = io_environment.defineParameterizationIO();

    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineClosure<WeaklyCompressibleFluid, ParameterizedViscosity>(
        ConstructArgs(rho0_f, c_f), ConstructArgs(parameterization_io, mu_f));
    water_block.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_locations);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicAlongAxis periodic_along_y(water_block.getSPHBodyBounds(), yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, periodic_along_x);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_block, periodic_along_y);

    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerNoRiemann> density_relaxation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationInner> update_density_by_summation(water_block_inner);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<NoLimiter, AllParticles>> transport_velocity_correction(water_block_inner);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    BodyRegionByCell free_stream_buffer(water_block, makeShared<MultiPolygonShape>(createBufferShape()));
    SimpleDynamics<FreeStreamCondition> freestream_condition(free_stream_buffer);


    NearShapeSurfaceStationaryBoundary near_surface_cylinder(water_block, makeShared<InverseShape<Cylinder>>("Cylinder"));
    near_surface_cylinder.getLevelSetShape().writeLevelSet(sph_system);
    fluid_dynamics::StationaryConfinementSimpleMethod confinement_condition_cylinder(near_surface_cylinder);

    update_density_by_summation.post_processes_.push_back(&confinement_condition_cylinder.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_cylinder.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_cylinder.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_cylinder.surface_bounding_);
    transport_velocity_correction.post_processes_.push_back(&confinement_condition_cylinder.transport_velocity_);
    viscous_force.post_processes_.push_back(&confinement_condition_cylinder.viscous_force_);
    
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);

    /** Computing viscous force acting on wall with wall model. */
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force_from_fluid(water_block, "ViscousForceFromFluid");
    //ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_from_fluid(cylinder, "PressureForceFromFluid");
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force_on_wall(water_block, "ViscousForceOnWall");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force(water_block, "ViscousForce");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_kernel_gradient_on_wall(water_block, "KernelGradientWall");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_kernel_gradient(water_block, "KernelGradient");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_force(water_block, "Force");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_force_inner(water_block, "ViscousForceInner");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 200.0;
    Real output_interval = end_time / 200.0;

    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;

        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;

            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                //pressure_force_on_cylinder.exec();
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                freestream_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            /** one need update configuration after periodic condition. */
            water_block_inner.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        write_total_viscous_force_on_wall.writeToFile(number_of_iterations);
        write_total_viscous_force.writeToFile(number_of_iterations);
        write_total_kernel_gradient_on_wall.writeToFile(number_of_iterations);
        write_total_kernel_gradient.writeToFile(number_of_iterations);
        write_total_force.writeToFile(number_of_iterations);
        write_total_force_inner.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    //write_total_viscous_force_on_inserted_body.testResult();

    return 0;
}
