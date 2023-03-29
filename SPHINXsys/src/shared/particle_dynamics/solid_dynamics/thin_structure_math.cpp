#include "thin_structure_math.h"

namespace SPH
{
	//=====================================================================================================//
	namespace thin_structure_dynamics
	{
		//=================================================================================================//
		Vec2d getVectorAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles)
		{
			/**The rotation matrix. */
			Real sin_angle = sin(rotation_angles[0]);
			Real cos_angle = cos(rotation_angles[0]);
			
			Mat2d rotation_matrix{
								  {cos_angle, sin_angle}, 	// First row
								  {-sin_angle,cos_angle},	//Second row
			};

			return rotation_matrix * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles)
		{
			/**The rotation matrix is the rotation around Y-axis \times rotation around X-axis. */
			Real sin_angle_x = sin(rotation_angles[0]);
			Real cos_angle_x = cos(rotation_angles[0]);
			
			Real sin_angle_y = sin(rotation_angles[1]);
			Real cos_angle_y = cos(rotation_angles[1]);

			Mat3d rotation_matrix = Mat3d{
				{cos_angle_y,  0.0, sin_angle_y},
				{0.0, 1.0, 0.0},
				{-sin_angle_y, 0.0, cos_angle_y},
			} * Mat3d{
				{1.0, 0.0, 0.0},
				{0.0, cos_angle_x, -sin_angle_x},
				{0.0, sin_angle_x,  cos_angle_x},
			};

			return rotation_matrix * initial_vector;
		}
		//=================================================================================================//
		Vec2d getVectorChangeRateAfterThinStructureRotation(const Vec2d &initial_vector, const Vec2d &rotation_angles, const Vec2d &angular_vel)
		{
			/**The derivative of the rotation matrix. */
			Real sin_angle = sin(rotation_angles[0]);
			Real cos_angle = cos(rotation_angles[0]);
			Mat2d drotation_matrix_dt{
				{-sin_angle * angular_vel[0],  cos_angle * angular_vel[0]},
				{-cos_angle * angular_vel[0], -sin_angle * angular_vel[0]},
			};

			return drotation_matrix_dt * initial_vector;
		}
		//=================================================================================================//
		Vec3d getVectorChangeRateAfterThinStructureRotation(const Vec3d &initial_vector, const Vec3d &rotation_angles, const Vec3d &angular_vel)
		{
			/**The rotation matrix about the X-axis. */
			Real sin_angle_x = sin(rotation_angles[0]);
			Real cos_angle_x = cos(rotation_angles[0]);
			
			Real sin_angle_y = sin(rotation_angles[1]);
			Real cos_angle_y = cos(rotation_angles[1]);

			Real angular_vel_x = angular_vel[0];
			Real angular_vel_y = angular_vel[1];

			Mat3d rotation_matrix_x{
				{1.0, 0.0, 0.0},
				{0.0, cos_angle_x, -sin_angle_x},
				{0.0, sin_angle_x, cos_angle_x},
			};
			/**The rotation matrix about the Y-axis. */
			Mat3d rotation_matrix_y{
				{cos_angle_y, 0.0, sin_angle_y},
				{0.0, 1.0, 0.0},
				{-sin_angle_y, 0.0, cos_angle_y},
			};
			/**The derivative of the rotation matrix of the X-axis. */
			Mat3d drotation_matrix_x_dt{
				{0.0, 0.0, 0.0},
				{0.0,-sin_angle_x * angular_vel_x, -cos_angle_x * angular_vel_x},
				{0.0, cos_angle_x * angular_vel_x, -sin_angle_x * angular_vel_x},
			};
			/**The derivative of the rotation matrix of the Y-axis. */
			Mat3d drotation_matrix_y_dt{
				{-sin_angle_y * angular_vel_y, 0.0, cos_angle_y * angular_vel_y},
				{0.0, 0.0, 0.0},
				{-cos_angle_y * angular_vel_y, 0.0,-sin_angle_y * angular_vel_y},
			};

			return (drotation_matrix_y_dt * rotation_matrix_x + rotation_matrix_y * drotation_matrix_x_dt)* initial_vector;
		}
		//=================================================================================================//
		Vec2d getRotationFromPseudoNormalForFiniteDeformation(const Vec2d &dpseudo_n_d2t, const Vec2d &rotation, const Vec2d &angular_vel, Real dt)
		{
			Vec2d dangular_vel_dt = Vec2d::Zero();
			dangular_vel_dt[0] = -(dpseudo_n_d2t[0] + sin(rotation[0]) * pow(angular_vel[0], 2))
								 / (2 * sin(rotation[0]) * angular_vel[0] * dt - cos(rotation[0]));
			return dangular_vel_dt;
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForFiniteDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
		{
			Real sin_rotation_x = sin(rotation[0]);
			Real cos_rotation_x = cos(rotation[0]);
			Real sin_rotation_y = sin(rotation[1]);
			Real cos_rotation_y = cos(rotation[1]);

			Real angle_vel_dt_x = (dpseudo_n_d2t[1] - sin_rotation_x * pow(angular_vel[0], 2))
								 / (2 * sin_rotation_x * angular_vel[0] * dt - cos_rotation_x);

			Real angle_vel_dt_y = (dpseudo_n_d2t[0] + cos_rotation_x * sin_rotation_y
								  * (pow(angular_vel[0], 2) + pow(angular_vel[1], 2))
								  + 2 * sin_rotation_x * cos_rotation_y * angular_vel[0] * angular_vel[1]
								  + (2 * cos_rotation_x * sin_rotation_y * angular_vel[0] * dt
								  + 2 * sin_rotation_x * cos_rotation_y * angular_vel[1] * dt
								  + sin_rotation_x * cos_rotation_y) * angle_vel_dt_x)
								 / (-2 * sin_rotation_x * cos_rotation_y * angular_vel[0] * dt
									- 2 * cos_rotation_x * sin_rotation_y * angular_vel[1] * dt
									+ cos_rotation_x * cos_rotation_y);
									
			return Vec3d(angle_vel_dt_x, angle_vel_dt_y, 0.0);
		}
		//=================================================================================================//
		Vec2d getRotationFromPseudoNormalForSmallDeformation(const Vec2d &dpseudo_n_d2t, const Vec2d &rotation, const Vec2d &angular_vel, Real dt)
		{
			return Vec2d(dpseudo_n_d2t[0], 0);
		}
		//=================================================================================================//
		Vec3d getRotationFromPseudoNormalForSmallDeformation(const Vec3d &dpseudo_n_d2t, const Vec3d &rotation, const Vec3d &angular_vel, Real dt)
		{
			return Vec3d(-dpseudo_n_d2t[1], dpseudo_n_d2t[0], 0.0);
		}
		//=================================================================================================//
		Vec2d getNormalFromDeformationGradientTensor(const Mat2d &F)
		{
			return Vec2d(-F.col(0)[1], F.col(0)[0]).normalized();
		}
		//=================================================================================================//
		Vec3d getNormalFromDeformationGradientTensor(const Mat3d &F)
		{
			return F.col(0).cross(F.col(1)).normalized();
		}
		//=================================================================================================//
		Vecd getLinearVariableJump(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			return particle_i_value - particle_j_value
				   - 0.5 * r_ij * (gradient_particle_i_value + gradient_particle_j_value) * e_ij;
		}
		//=================================================================================================//
		Vecd getWENOVariableJump(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			return getWENOLeftState(e_ij, r_ij, particle_i_value,
									gradient_particle_i_value, particle_j_value, gradient_particle_j_value)
				   - getWENORightState(e_ij, r_ij, particle_i_value,
									   gradient_particle_i_value, particle_j_value, gradient_particle_j_value);
		}
		//=================================================================================================//
		Vecd getWENOStateWithStencilPoints(const Vecd &v1, const Vecd &v2, const Vecd &v3, const Vecd &v4)
		{
			Vecd f1 = 0.5 * v2 + 0.5 * v3;
			Vecd f2 = -0.5 * v1 + 1.5 * v2;
			Vecd f3 = v2 / 3.0 + 5.0 * v3 / 6.0 - v4 / 6.0;

			Real epsilon = 1.0e-6;
			Real s1 = (v2 - v3).dot(v2 - v3) + epsilon;
			Real s2 = (v2 - v1).dot(v2 - v1) + epsilon;
			Real s3 = (3.0 * v2 - 4.0 * v3 + v4).dot(3.0 * v2 - 4.0 * v3 + v4) / 4.0
					  + 13.0 * (v2 - 2.0 * v3 + v4).dot(v2 - 2.0 * v3 + v4) / 12.0 + epsilon;
			Real s12 = 13.0 * (v1 - 2.0 * v2 + v3).dot(v1 - 2.0 * v2 + v3) / 12.0
					   + (v1 - v3).dot(v1 - v3) / 4.0 + epsilon;
			Real s4 = (v1.dot(6649.0 * v1 - 30414.0 * v2 + 23094.0 * v3 - 5978.0 * v4)
					   + 3.0 * v2.dot(13667.0 * v2 - 23534.0 * v3 + 6338.0 * v4)
					   + 3.0 * v3.dot(11147.0 * v3 - 6458.0 * v4)
					   + 3169.0 * v4.dot(v4)) / 2880.0;
			Real tau_4 = s4 - 0.5 * (s1 + s2);

			Real alpha_1 = (1.0 + (tau_4 / s1) * (tau_4 / s12)) / 3.0;
			Real alpha_2 = (1.0 + (tau_4 / s2) * (tau_4 / s12)) / 6.0;
			Real alpha_3 = (1.0 + tau_4 / s3) / 2.0;
			Real w_1 = alpha_1 / (alpha_1 + alpha_2 + alpha_3);
			Real w_2 = alpha_2 / (alpha_1 + alpha_2 + alpha_3);
			Real w_3 = alpha_3 / (alpha_1 + alpha_2 + alpha_3);

			return w_1 * f1 + w_2 * f2 + w_3 * f3;
		}
		//=================================================================================================//
		Vecd getWENOLeftState(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			Vecd v1 = particle_i_value + gradient_particle_i_value * e_ij * r_ij;
			Vecd v2 = particle_i_value;
			Vecd v3 = particle_j_value;
			Vecd v4 = particle_j_value - gradient_particle_j_value * e_ij * r_ij;

			return getWENOStateWithStencilPoints(v1, v2, v3, v4);
		}
		//=================================================================================================//
		Vecd getWENORightState(const Vecd &e_ij, const Real &r_ij, const Vecd &particle_i_value,
			const Matd &gradient_particle_i_value, const Vecd &particle_j_value, const Matd &gradient_particle_j_value)
		{
			Vecd v1 = particle_j_value - gradient_particle_j_value * e_ij * r_ij;
			Vecd v2 = particle_j_value;
			Vecd v3 = particle_i_value;
			Vecd v4 = particle_i_value + gradient_particle_i_value * e_ij * r_ij;

			return getWENOStateWithStencilPoints(v1, v2, v3, v4);
		}
		//=================================================================================================//
		Vec2d getRotationJump(const Vec2d &pseudo_n_jump, const Mat2d &transformation_matrix)
		{
			Vec2d local_rotation_jump =  Vec2d::Zero();
			Vec2d local_pseuodo_n_jump = transformation_matrix * pseudo_n_jump;
			local_rotation_jump[0] = local_pseuodo_n_jump[0];
			return transformation_matrix.transpose() * local_rotation_jump;
		}
		//=================================================================================================//
		Vec3d getRotationJump(const Vec3d &pseudo_n_jump, const Mat3d &transformation_matrix)
		{
			Vec3d local_rotation_jump = Vec3d::Zero();
			Vec3d local_pseuodo_n_jump = transformation_matrix * pseudo_n_jump;
			local_rotation_jump[0] = local_pseuodo_n_jump[0];
			local_rotation_jump[1] = local_pseuodo_n_jump[1];
			return transformation_matrix.transpose() * local_rotation_jump;
		}
		//=================================================================================================//
		Mat2d getCorrectedAlmansiStrain(const Mat2d &current_local_almansi_strain, const Real &nu_)
		{
			Mat2d corrected_almansi_strain = current_local_almansi_strain;
			corrected_almansi_strain(1,1) = -nu_ * current_local_almansi_strain(0,0) / (1.0 - nu_);
			return corrected_almansi_strain;
		}
		//=================================================================================================//
		Mat3d getCorrectedAlmansiStrain(const Mat3d &current_local_almansi_strain, const Real &nu_)
		{
			Mat3d corrected_almansi_strain = current_local_almansi_strain;
			corrected_almansi_strain(2,2)
				= -nu_ * (current_local_almansi_strain(0,0) + current_local_almansi_strain(1,1)) / (1.0 - nu_);
			return corrected_almansi_strain;
		}
		//=================================================================================================//
		Mat2d getCorrectionMatrix(const Mat2d &local_deformation_part_one)
		{
			Real one_over_local_deformation = 1.0 / local_deformation_part_one(0,0);
			return Mat2d{
							{one_over_local_deformation, 0},
							{0,0},
						};
		}
		//=================================================================================================//
		Mat3d getCorrectionMatrix(const Mat3d &local_deformation_part_one)
		{
			Mat3d correction_matrix = Mat3d::Zero();
			correction_matrix.block<2,2>(0,0) = local_deformation_part_one.block<2,2>(0,0).inverse();
			return correction_matrix;
		}
		//=================================================================================================//
	}
}
