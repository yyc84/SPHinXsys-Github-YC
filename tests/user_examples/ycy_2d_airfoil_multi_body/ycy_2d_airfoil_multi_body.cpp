/**
* @file 	particle_generator_single_resolution.cpp
* @brief 	This is the test of using level set to generate particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by level set for a complex geometry (2D).
*			Before particle generation, we clean the level set, then do re-initialization.

* @author 	Yongchuan Yu and Xiangyu Hu
*/

#include "sphinxsys.h"

using namespace SPH;

std::string airfoil_flap_front = "./input/airfoil_flap_front.dat";
std::string airfoil_wing = "./input/airfoil_wing.dat";
std::string airfoil_flap_rear = "./input/airfoil_flap_rear.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.25;             /**< airfoil length rear part. */
Real DL1 = 0.25;            /**< airfoil length front part. */
Real DH = 0.25;             /**< airfoil height. */
Real resolution_ref = 0.0025; /**< Reference resolution. */
BoundingBox system_domain_bounds(Vec2d(-0.3*DL, -1.5*DH), Vec2d(1.0*DL, 1.5 * DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * (DL + DL1) ) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Shape of the InputBody
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-0.3*DL, -1.5*DH));
    water_block_shape.push_back(Vecd(-0.3*DL, 1.5* DH ));
    water_block_shape.push_back(Vecd(1.0*DL, 1.5* DH));
    water_block_shape.push_back(Vecd(1.0*DL, -1.5*DH));
    water_block_shape.push_back(Vecd(-0.3*DL, -1.5*DH));

    return water_block_shape;
}

class ImportModel : public MultiPolygonShape
{
  public:
    explicit ImportModel(const std::string &import_model_name) : MultiPolygonShape(import_model_name)
    {
        multi_polygon_.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
    }
};


class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon airfoil;
        airfoil.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
        airfoil.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
        airfoil.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
        subtract<MultiPolygonShape>(airfoil);
    }
};

//class WaterShape : public ComplexLevelSetShape
//{
//  public:
//    explicit WaterShape(const std::string &shape_name) : ComplexLevelSetShape(shape_name)
//    {
//        MultiPolygon outer_boundary(createWaterBlockShape());
//        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
//        MultiPolygon airfoil;
//        airfoil.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
//        airfoil.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
//        airfoil.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
//        subtract<MultiPolygonShape>(airfoil);
//    }
//};
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    BOOLEAN complex_relaxation(true);
    //SPHAdaptation sph_adaptation(resolution_ref);
    //WaterShape water_shape("WaterShape");
    /*water_shape.initializeComponentLevelSetShapesByAdaptation(makeShared<SPHAdaptation>(resolution_ref), io_environment);
    StdVec<LevelSetShape*> shapes(water_shape.getLevelSetShapes());
    for (size_t i = 0; i != water_shape.getLevelSetShapes().size(); ++i)
    {
        water_shape.getLevelSetShapes()[i]->writeLevelSet(io_environment);
    }*/
    
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    //RealBody input_body(sph_system, makeShared<InputBody>("SPHInXsysLogo"));
    RealBody airfoil(sph_system, makeShared<ImportModel>("AirFoil"));
    airfoil.defineBodyLevelSetShape()->cleanLevelSet()->writeLevelSet(io_environment);
    //airfoil.initialComponentLevelSetShapes();
    airfoil.defineParticlesAndMaterial();
    airfoil.generateParticles<ParticleGeneratorLattice>();

 

    RealBody water(sph_system, makeShared<WaterBlock>("Water"));
    water.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(io_environment);
    //water.defineComponentLevelSetShape("AirFoil")->cleanLevelSet()->writeLevelSet(io_environment);
    water.defineParticlesAndMaterial();
    water.generateParticles<ParticleGeneratorLattice>();


    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation airfoil_inner(airfoil);
    InnerRelation water_inner(water);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
   
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt airfoil_recording_to_vtp(io_environment, {&airfoil});
    MeshRecordingToPlt cell_linked_list_recording(io_environment, airfoil.getCellLinkedList());

    BodyStatesRecordingToPlt water_recording_to_vtp(io_environment, { &water });
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    if (complex_relaxation == 0)
    {
      using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_airfoil_particles(airfoil);
        relax_dynamics::RelaxationStepInner relaxation_step_inner_airfoil(airfoil_inner, true);

        SimpleDynamics<RandomizeParticlePosition> random_water_particles(water);
        relax_dynamics::RelaxationStepInner relaxation_step_inner_water(water_inner, true);
        random_airfoil_particles.exec(0.25);
        relaxation_step_inner_airfoil.SurfaceBounding().exec();
        airfoil.updateCellLinkedList();
        random_water_particles.exec(0.25);
        relaxation_step_inner_water.SurfaceBounding().exec();
        water.updateCellLinkedList();
        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        airfoil_recording_to_vtp.writeToFile();
        cell_linked_list_recording.writeToFile();
        water_recording_to_vtp.writeToFile();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner_airfoil.exec();
            relaxation_step_inner_water.exec();
            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                airfoil_recording_to_vtp.writeToFile(ite_p);
                water_recording_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;

        return 0;
    }
    
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    //water_block.sph_adaptation_->resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    water_block.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(io_environment);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    BodyStatesRecordingToVtp water_block_recording_to_vtp(io_environment, { &water_block });

    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&airfoil});
    ContactRelation airfoil_contact(airfoil, {&water_block});
    ComplexRelation water_airfoil_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    //InnerRelation cylinder_inner(cylinder);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    using namespace relax_dynamics;
    SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(airfoil);
    SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
    BodyStatesRecordingToPlt write_real_body_states(io_environment, sph_system.real_bodies_);
    ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);
    relax_dynamics::RelaxationStepInner relaxation_step_inner(airfoil_inner);
    relax_dynamics::RelaxationStepComplex relaxation_step_complex(water_airfoil_complex, "OuterBoundary", true);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_inserted_body_particles.exec(0.25);
    random_water_body_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    relaxation_step_complex.SurfaceBounding().exec();
    write_real_body_states.writeToFile(0);
    //airfoil_recording_to_vtp.writeToFile();
    //water_recording_to_vtp.writeToFile();
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        relaxation_step_complex.exec();
        ite_p += 1;
        if (ite_p % 200 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            write_real_body_states.writeToFile(ite_p);
            //airfoil_recording_to_vtp.writeToFile(ite_p);
            //water_recording_to_vtp.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;

    //write_real_body_particle_reload_files.writeToFile(0);

    return 0;
}
