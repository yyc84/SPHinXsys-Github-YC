
#include "chrono/physics/ChSystemSMC.h"

#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "chrono/solver/ChSolverPMINRES.h" 
#include "chrono/solver/ChSolverPSOR.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>

using namespace chrono;
using namespace chrono::geometry;


void WriteMeshVTP(const std::shared_ptr<ChTriangleMeshConnected>& mesh,
                  const ChFrame<>& X, const std::string& fname) {
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
           "<Piece NumberOfPoints=\"" << Nv << "\" NumberOfPolys=\"" << Nf << "\">\n"
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

    if (!V.empty()) {
        auto v_local = V[0];
        auto v_world = X.TransformPointLocalToParent(v_local);
        std::cout << "[DEBUG VTP] First vertex local: " << v_local
            << ", world: " << v_world
            << ", frame pos: " << X.GetPos() << std::endl;
    }
}

int main() {
    std::filesystem::path cwd = std::filesystem::current_path();
    std::filesystem::path input_dir = cwd / "input";
    
    std::filesystem::path cube_file = input_dir / "cube03.obj";
    std::filesystem::path sphere_file = input_dir / "sphere03.obj";
    std::string cubePath = cube_file.string();
    std::string spherePath = sphere_file.string();

    std::filesystem::path outdir = std::filesystem::current_path() / "output";
    std::filesystem::create_directories(outdir);

    ChSystemSMC system;
    auto solver = chrono_types::make_shared<ChSolverPSOR>();
    solver->SetMaxIterations(100);   
    //solver->SetTolerance(1e-10);    
    //solver->EnableDiagonalPreconditioner(true); 
    system.SetSolver(solver);  
    // 其他设置
    system.Set_G_acc(ChVector<>(0, -9.81, 0));
    system.SetMaxPenetrationRecoverySpeed(1.0);

    std::cout << "Solver type: " << typeid(*solver).name() << std::endl;

    std::ifstream test(cubePath);
    if (!test.good()) {
        std::cerr << "ERROR: File not found at: " << cubePath << std::endl;
        return 1;
    }

    auto cube_mesh = chrono_types::make_shared<chrono::geometry::ChTriangleMeshConnected>();
    if (!cube_mesh->LoadWavefrontMesh(cubePath, false)) {
        std::cerr << "Failed to load cube.stl!" << std::endl;
        return 1;
    } else {
        const auto& V = cube_mesh->getCoordsVertices();    
        const auto& F = cube_mesh->getIndicesVertexes();
        std::cout << "Loaded mesh: vertices=" << V.size()
                  << ", triangles=" << F.size() << "\n ";
    }

    auto cube_body = chrono_types::make_shared<ChBody>();
    cube_body->SetBodyFixed(true);
    cube_body->SetCollide(true);

    auto cube_material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    cube_material->SetFriction(0.4f);      
    cube_material->SetRestitution(0.1f);   
    cube_material->SetYoungModulus(1e7f);  
    cube_material->SetAdhesion(0);
    cube_material->SetKn(1e6f);
    cube_material->SetKt(1e5f);
    cube_material->SetGn(1e2f);
    cube_material->SetGt(1e2f);

    cube_body->GetCollisionModel()->ClearModel();
    cube_body->GetCollisionModel()->AddTriangleMesh(
        cube_material,
        cube_mesh, false, false,
        ChVector<>(0, 0, 0),
        ChMatrix33<>(1) 
    );
    cube_body->GetCollisionModel()->BuildModel();
    cube_body->SetCollide(true);

    auto mesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
    mesh_shape->SetMesh(cube_mesh);
    mesh_shape->SetBackfaceCull(true);  // 可选：启用背面剔除
    cube_body->AddVisualShape(mesh_shape);
    system.Add(cube_body);

    auto sphere_mesh = chrono_types::make_shared<chrono::geometry::ChTriangleMeshConnected>();
    //sphere_mesh->LoadSTLMesh(spherePath, false);
    if (!sphere_mesh->LoadWavefrontMesh(spherePath, false)) {
        std::cerr << "Failed to load cube.stl!" << std::endl;
        return 1;
    } else {
        const auto& V = sphere_mesh->getCoordsVertices();
        const auto& F = sphere_mesh->getIndicesVertexes();
        std::cout << "Loaded mesh: vertices=" << V.size()
                  << ", triangles=" << F.size() << "\n ";
    }

    auto sphere_body = chrono_types::make_shared<ChBody>();
    sphere_body->SetMass(1.0);
    sphere_body->SetInertiaXX(ChVector<>(1, 1, 1));
    sphere_body->SetPos(ChVector<>(0, 0, 0.0));
    sphere_body->SetCollide(true);

    auto sphere_material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    sphere_material->SetFriction(0.4f);
    sphere_material->SetRestitution(0.1f);
    sphere_material->SetYoungModulus(1e7f);
    sphere_material->SetAdhesion(0);
    sphere_material->SetKn(1e6f);
    sphere_material->SetKt(1e5f);
    sphere_material->SetGn(1e2f);
    sphere_material->SetGt(1e2f);

    sphere_body->GetCollisionModel()->ClearModel();
    sphere_body->GetCollisionModel()->AddTriangleMesh(
        sphere_material,
        sphere_mesh,
        false,                  // is_static
        false,                  // is_convex
        ChVector<>(0, 0, 0),
        ChMatrix33<>(1)
    );
    sphere_body->GetCollisionModel()->BuildModel();
    sphere_body->SetCollide(true);

    auto sphere_visual = chrono_types::make_shared<ChTriangleMeshShape>();
    sphere_visual->SetMesh(sphere_mesh);
    sphere_visual->SetBackfaceCull(true);
    sphere_body->AddVisualShape(sphere_visual);

    system.Add(sphere_body);

    //std::filesystem::create_directory("output");
    WriteMeshVTP(cube_mesh, ChFrame<>(), (outdir / "cube_initial.vtp").string());
    WriteMeshVTP(sphere_mesh, ChFrame<>(), (outdir / "sphere_initial.vtp").string());

    double step_size = 1e-4;
    double t_end = 0.1;
    int step = 0;
    int out_every = 5;

    while (system.GetChTime() < t_end) {
        system.DoStepDynamics(step_size);

        if (step % out_every == 0) {
            std::cout << "[t=" << system.GetChTime() << "] Sphere z = " << sphere_body->GetPos().z() << std::endl;

            std::string base = "output/mesh_" + std::to_string(step);
            WriteMeshVTP(cube_mesh, ChFrame<>(cube_body->GetPos(), cube_body->GetRot()), base + "_cube.vtp");
            WriteMeshVTP(sphere_mesh, ChFrame<>(sphere_body->GetPos(), sphere_body->GetRot()), base + "_sphere.vtp");

            if (sphere_body->GetPos().z() < 0.3 && sphere_body->GetPos_dt().z() < 0.0) {
                std::cout << "[t=" << system.GetChTime() << "] Contact likely happening.\n";
            }

            std::cout << "[DEBUG] sphere_body pos z = " << sphere_body->GetPos().z() << std::endl;

        }

        ++step;
    }

    return 0;
}
