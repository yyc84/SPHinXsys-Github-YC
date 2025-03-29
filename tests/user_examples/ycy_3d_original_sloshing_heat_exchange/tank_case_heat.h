/*
* @file 	waterentry_elastic_case.h
*/

#ifndef	TANK_CASE_HEAT_H
#define TANK_CASE_HEAT_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/*@brief Basic geometry parameters and numerical setup.
*/


Real resolution_ref = 0.008;   /* Initial particle spacing*/
/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.2, -0.05, -0.2), Vec3d(0.2, 1.0, 0.2));


/*Material properties of the water and air.*/
//Real rho0_f = 1000.0;         /*Fluid density*/
//Real rho0_a = 1.226;            /*Air density*/
//Real gravity_g = 9.81;        /*Gravity force of fluid*/
//Real U_f = 2.0 * sqrt(gravity_g * 0.5);	/**< Characteristic velocity. */
//Real U_g = 2.0 * sqrt(gravity_g * 0.5);  	/**< dispersion velocity in shallow water. */
//Real c_f = 10.0 * SMAX(U_g, U_f);	/**< Reference sound speed. */
//Real f = 1.7;
//Real a = 0.04;
//Real c_p_water = 4.179e3;
//Real c_p_air = 1.012e3;
//Real k_water = 0.620;
//Real k_air = 0.0254;
//Real diffusion_coff_water = k_water / (c_p_water * rho0_f);
//Real diffusion_coff_air = k_air / (c_p_air * rho0_a);
//Real mu_water = 653.9e-6;
//Real mu_air = 20.88e-6;
//Real length_scale = 1.0;
//Vec3d translation(0, 0.0, 0);
//Real initial_temperature_right = 313.15;
//Real initial_temperature_air = 353.15;

/*parameters for Methane and Nitrogen*/
Real rho0_f = 424.7;                   /*Fluid density*/
Real rho0_a = 2.668;                    /*Air density*/
Real gravity_g = 9.81;                  /*Gravity force of fluid*/
Real U_f = 2.0 * sqrt(gravity_g * 0.5); /**< Characteristic velocity. */
Real U_g = 2.0 * sqrt(gravity_g * 0.5); /**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);       /**< Reference sound speed. */
Real f = 1.3;
Real a = 0.02;
Real c_p_water = 3.4267e3;
Real c_p_air = 1.054e3;
Real k_water = 0.1846;
Real k_air = 0.01221;
Real diffusion_coff_water = k_water / (c_p_water * rho0_f);
Real diffusion_coff_air = k_air / (c_p_air * rho0_a);
Real mu_water = 121.79e-6;
Real mu_air = 8.73e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.0, 0);
Real initial_temperature_right = 110.15;
Real initial_temperature_air = 130.15;


std::string fuel_tank_outer = "./input/tank_outer_small.stl";
std::string air_05 = "./input/gas_small.stl";
std::string probe_s1_shape = "./input/ProbeS1.stl";
std::string probe_s2_shape = "./input/ProbeS2.stl";
std::string probe_s3_shape = "./input/ProbeS3.stl";
//std::string fuel_tank_inner = "./input/tank_inner_small.stl";
//std::string water_05 = "./input/water_small.stl";

//std::string fuel_tank_inner = "./input/tank_inner_2hole_small.stl";
//std::string water_05 = "./input/water_2hole_small.stl";

std::string fuel_tank_inner = "./input/tank_inner_9hole_small.stl";
std::string water_05 = "./input/water_9hole_small.stl";

//std::string fuel_tank_inner = "./input/tank_inner_ring_small.stl";
//std::string water_05 = "./input/water_ring_small.stl";

/*Fuel Tank.*/
class Tank : public ComplexShape
{
public:
	explicit Tank(const std::string& shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */

		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale, "OuterWall");
		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale, "InnerWall");
	}
};
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};

class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
	}
};

class VariableGravity : public Gravity
{
    Real time_ = 0;

  public:
    VariableGravity(Vecd gravity_vector) : Gravity(gravity_vector){};
    virtual Vecd InducedAcceleration(const Vecd &position = Vecd::Zero()) override
    {
        time_ = GlobalStaticVariables::physical_time_;
        if (time_ <= 2.0)
        {
            global_acceleration_ = global_acceleration_;
        }
        else
        {
            global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * (time_-2));
        }
        //global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
        return global_acceleration_;
    }
};

class ProbeS1 : public ComplexShape
{
public:
	explicit ProbeS1(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
	}
};

class ProbeS2 : public ComplexShape
{
public:
	explicit ProbeS2(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_2(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
	}
};

class ProbeS3 : public ComplexShape
{
public:
	explicit ProbeS3(const std::string& shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_3(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
	}
};

class ThermoAirBodyInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermoAirBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")), 
          heat_flux_inner_(*particles_->registerSharedVariable<Real>("PhiFluxInner")),
          heat_flux_contact_(*particles_->registerSharedVariable<Real>("PhiFluxContact")),
          heat_flux_wu_(*particles_->registerSharedVariable<Real>("PhiFluxWuContact")) {};
    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_air;
    };

  protected:
    StdLargeVec<Real> &phi_;
    StdLargeVec<Real> &heat_flux_inner_;
    StdLargeVec<Real> &heat_flux_contact_;
    StdLargeVec<Real> &heat_flux_wu_;
};

class ThermoWaterBodyInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit ThermoWaterBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
          heat_flux_inner_(*particles_->registerSharedVariable<Real>("PhiFluxInner")),
          heat_flux_contact_(*particles_->registerSharedVariable<Real>("PhiFluxContact")), 
          heat_flux_wu_(*particles_->registerSharedVariable<Real>("PhiFluxWuContact")) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature_right;
    };

  protected:
    StdLargeVec<Real> &phi_;
    StdLargeVec<Real> &heat_flux_inner_;
    StdLargeVec<Real> &heat_flux_contact_;
    StdLargeVec<Real> &heat_flux_wu_;
};

using HeatExchangeComplex = HeatExchangeDiffusionComplex<KernelGradientInner, KernelGradientContact, HeatIsotropicDiffusion, HeatIsotropicDiffusion>;

StdVec<Vecd> LiquidTemperatureObserverParticle()
{
    StdVec<Vecd> observation_points;

    observation_points.push_back(Vecd(0.0, 0.35, 0.0));
    observation_points.push_back(Vecd(0.0, 0.5, 0.0));
    observation_points.push_back(Vecd(0.0, 0.6, 0.0));

    return observation_points;
};

StdVec<Vecd> GasTemperatureObserverParticle()
{
    StdVec<Vecd> observation_points;

	observation_points.push_back(Vecd(0.0, 0.4, 0.0));
    observation_points.push_back(Vecd(0.0, 0.5, 0.0));
    observation_points.push_back(Vecd(0.0, 0.65, 0.0));

    return observation_points;
};

#endif //TANK_CASE_HEAT_H