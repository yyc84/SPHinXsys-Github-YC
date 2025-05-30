#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/solver/ChSolverPSOR.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/physics/ChInertiaUtils.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>


using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;
using namespace chrono::collision;

void WriteMeshVTP(const std::shared_ptr<ChTriangleMeshConnected>& mesh, const ChFrame<>& X, const std::string& fname) {
    std::ofstream out(fname);
    if (!out.good()) {
        std::cerr << "Cannot write to " << fname << std::endl;
        return;
    }

    const auto& V = mesh->getCoordsVertices();
    const auto& F = mesh->getIndicesVertexes();

    int Nv = static_cast<int>(V.size());
    int Nf = static_cast<int>(F.size());

    out << "<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           "<PolyData>\n"
           "<Piece NumberOfPoints=\""
        << Nv << "\" NumberOfPolys=\"" << Nf
        << "\">\n"
           "<Points>\n"
           "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& v : V) {
        auto p = X.TransformPointLocalToParent(v);
        out << "  " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    out << "</DataArray>\n</Points>\n"
           "<Polys>\n"
           "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for (const auto& f : F) {
        out << "  " << f.x() << " " << f.y() << " " << f.z() << "\n";
    }

    out << "</DataArray>\n"
           "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

    for (int i = 1; i <= Nf; ++i)
        out << "  " << i * 3 << "\n";

    out << "</DataArray>\n</Polys>\n"
           "</Piece>\n"
           "</PolyData>\n"
           "</VTKFile>\n";

    out.close();
}

class MyContactReporter : public chrono::ChContactContainer::ReportContactCallback {
  public:
    virtual bool OnReportContact(const ChVector<>& pA,             // 接触点 A
                                 const ChVector<>& pB,             // 接触点 B
                                 const ChMatrix33<>& plane_coord,  // 接触坐标系，Z轴是法线
                                 const double& distance,           // 分离距离（负数表示穿透）
                                 const double& eff_radius,         // 有效曲率半径
                                 const ChVector<>& react_forces,   // 反力
                                 const ChVector<>& react_torques,  // 反扭矩
                                 ChContactable* contactobjA,       // 接触体 A
                                 ChContactable* contactobjB        // 接触体 B
                                 ) override {
        std::cout << "[Contact] pA: " << pA << ", pB: " << pB << ", dist: " << distance << ", force: " << react_forces
                  << std::endl;
        return true;  // 返回 true 表示继续报告下一个接触
    }
};
std::string cubePath = "./input/cube01.obj";
std::string spherePath = "./input/sphere01.obj";

int main() 
{

    std::cout << "Chrono version: " << CHRONO_VERSION << std::endl;
    ChSystemNSC system;
    system.Set_G_acc(ChVector<>(0, -9.81, 0));
    system.SetMaxPenetrationRecoverySpeed(1.0);
    system.SetCollisionSystemType(chrono::collision::ChCollisionSystemType::BULLET);
    system.SetMinBounceSpeed(0.1); 

    collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
    collision::ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    auto solver = chrono_types::make_shared<ChSolverPSOR>();
    solver->SetMaxIterations(200);
    solver->EnableWarmStart(true);
    system.SetSolver(solver);

    // Load cube (static ground)

    auto cube_mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(cubePath, false, true);
    cube_mesh->RepairDuplicateVertexes(1e-9);  
    cube_mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(1));

    // compute mass inertia from mesh
    double cube_mass;
    ChVector<> cube_cog;
    ChMatrix33<> cube_inertia;
    double cube_density = 1000;
    cube_mesh->ComputeMassProperties(true, cube_mass, cube_cog, cube_inertia);
    ChMatrix33<> cube_principal_inertia_rot;
    ChVector<> cube_principal_I;
    ChInertiaUtils::PrincipalInertia(cube_inertia, cube_principal_I, cube_principal_inertia_rot);

    // Create a shared visual model containing a visualizatoin mesh
    auto cube_mesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    cube_mesh_shape->SetMesh(cube_mesh);
    cube_mesh_shape->SetMutable(false);
    cube_mesh_shape->SetColor(ChColor(1.0f, 0.5f, 0.5f));
    cube_mesh_shape->SetBackfaceCull(true);

    /*auto cube_vis_model = chrono_types::make_shared<ChVisualModel>();
    cube_vis_model->AddShape(cube_mesh_shape);*/

    auto cube_body = chrono_types::make_shared<ChBody>();
    //cube_body->SetMass(cube_mass * cube_density);
    cube_body->SetBodyFixed(true);
    cube_body->SetCollide(true);

    auto cube_material = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    cube_body->GetCollisionModel()->ClearModel();
    //cube_body->GetCollisionModel()->AddTriangleMesh(cube_material, cube_mesh, false, false, VNULL, ChMatrix33<>(1), 0.002);
    cube_body->GetCollisionModel()->AddTriangleMesh(cube_material, cube_mesh, false, false);
    cube_body->GetCollisionModel()->BuildModel();

    cube_body->AddVisualShape(cube_mesh_shape);
    system.Add(cube_body);

    // Load sphere (falling object)
    auto sphere_mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(spherePath, false, true);
    sphere_mesh->RepairDuplicateVertexes(1e-9);  // 修复重顶点
    sphere_mesh->Transform(ChVector<>(0, 0.1, 0), ChMatrix33<>(1));

    // compute mass inertia from mesh
    double sphere_mass;
    ChVector<> sphere_cog;
    ChMatrix33<> sphere_inertia;
    double sphere_density = 1000;
    sphere_mesh->ComputeMassProperties(true, sphere_mass, sphere_cog, sphere_inertia);
    ChMatrix33<> sphere_principal_inertia_rot;
    ChVector<> sphere_principal_I;
    ChInertiaUtils::PrincipalInertia(sphere_inertia, sphere_principal_I, sphere_principal_inertia_rot);

    auto sphere_mesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    sphere_mesh_shape->SetMesh(sphere_mesh);
    sphere_mesh_shape->SetMutable(false);
    sphere_mesh_shape->SetColor(ChColor(1.0f, 0.5f, 1.0f));
    sphere_mesh_shape->SetBackfaceCull(true);

    auto sphere_body = chrono_types::make_shared<ChBody>();
    sphere_body->SetMass(sphere_mass * sphere_density);
    sphere_body->SetInertiaXX(sphere_density * sphere_principal_I);
    sphere_body->SetPos(ChVector<>(0, 0.0, 0));
    sphere_body->SetCollide(true);

    auto sphere_material = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    sphere_material->SetFriction(0.2f);

    sphere_body->GetCollisionModel()->ClearModel();
    sphere_body->GetCollisionModel()->AddTriangleMesh(sphere_material, sphere_mesh, false, false);
    sphere_body->GetCollisionModel()->BuildModel();

    sphere_body->AddVisualShape(sphere_mesh_shape);
    system.Add(sphere_body);

    // Output initial meshes
    WriteMeshVTP(cube_mesh, ChFrame<>(cube_body->GetPos(), cube_body->GetRot()),
                 (outdir / "cube_initial.vtp").string());
    WriteMeshVTP(sphere_mesh, ChFrame<>(sphere_body->GetPos(), sphere_body->GetRot()),
                 (outdir / "sphere_initial.vtp").string());

    // Create Irrlicht visualization system
    auto vis = chrono_types::make_shared<chrono::irrlicht::ChVisualSystemIrrlicht>();
    vis->AttachSystem(&system);
    vis->SetWindowSize(1024, 768);
    vis->SetWindowTitle("Sphere falling on mesh box");
    vis->Initialize();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0.5, 0.5, -1.5));
    vis->AddTypicalLights();
    vis->AddLightWithShadow(ChVector<>(1.5, 5.5, -2.5), ChVector<>(0, 0, 0), 3, 2.2, 7.2, 40, 512,
                            ChColor(0.8f, 0.8f, 1.0f));

  

    double step_size = 1e-5;
    double t_end = 0.3;
    int step = 0;
    int out_every = 5;

    // Simulation loop with visualization
    while (system.GetChTime() < t_end)
    //while (vis->Run()) 
    {
        //vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        //vis->Render();
        //vis->EndScene();

        system.DoStepDynamics(step_size);

        if (system.GetChTime() > t_end)
            break;

        if (step % out_every == 0) {
            //std::string base = "output/mesh_" + std::to_string(step);
            WriteMeshVTP(cube_mesh, ChFrame<>(cube_body->GetPos(), cube_body->GetRot()), ("output/mesh_cube_"  + std::to_string(step) +".vtp"));
            WriteMeshVTP(sphere_mesh, ChFrame<>(sphere_body->GetPos(), sphere_body->GetRot()), ("output/mesh_sphere_" + std::to_string(step) + ".vtp"));

            std::cout << "[t=" << system.GetChTime() << "] Sphere y = " << sphere_body->GetPos().y() << std::endl;
            std::cout << "Number of bodies: " << system.Get_bodylist().size() << std::endl;
            //auto contact_reporter = std::make_shared<MyContactReporter>();
            //system.GetContactContainer()->ReportAllContacts(contact_reporter);
        }
        ++step;
    }

    return 0;
}
