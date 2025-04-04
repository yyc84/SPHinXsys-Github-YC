/*
* @file 	waterentry_elastic_case.h
*/

#ifndef	TANK_CASE_H
#define TANK_CASE_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/*@brief Basic geometry parameters and numerical setup.
*/

Real resolution_ref = 0.0035;   /* Initial particle spacing*/


/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3,-0.3), Vec3d(0.3, 0.5,0.3));


/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;         /*Fluid density*/
Real rho0_a = 1.0;            /*Air density*/
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_f = 2.0* sqrt(gravity_g * 0.174);	/**< Characteristic velocity. */
Real U_g = 2.0* sqrt(gravity_g * 0.174);  	/**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);	/**< Reference sound speed. */
Real f = 1.63;
Real a = 0.0075;
Real mu_water = 653.9e-6;
Real mu_air = 20.88e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.175, 0);

std::string fuel_tank_outer = "./input/validation_tank_outer_slim.stl";
std::string fuel_tank_inner = "./input/validation_tank_inner.stl";
std::string water_05 = "./input/validation_water.stl";
std::string air_05 = "./input/validation_air.stl";
std::string probe_s1_shape = "./input/ProbeS1.stl";
std::string probe_s2_shape = "./input/ProbeS2.stl";
std::string probe_s3_shape = "./input/ProbeS3.stl";

/*
Fuel Tank.
*/
class Tank : public ComplexShape
{
public:
	explicit Tank(const std::string &shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */

		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale,"OuterWall");
		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale,"InnerWall");
	}
};
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};

class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
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
		time_= GlobalStaticVariables::physical_time_;
            if (time_ <= 3.0)
            {
                global_acceleration_ = global_acceleration_;
            }
            else
            {
                global_acceleration_[0] = -4.0 * PI * PI * f * f * a * sin(2 * PI * f * (time_ - 3));
            }
            // global_acceleration_[0] = 4.0 * PI * PI * f * f * a * sin(2 * PI * f * time_);
            return global_acceleration_;
	}
};

class ProbeS1 : public ComplexShape
{
public:
	explicit ProbeS1(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
	}
};

class ProbeS2 : public ComplexShape
{
public:
	explicit ProbeS2(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_2(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
	}
};

class ProbeS3 : public ComplexShape
{
public:
	explicit ProbeS3(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_3(0.0,0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
	}
};

#endif //TANK_CASE_H