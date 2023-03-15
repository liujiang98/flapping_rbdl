/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <cmath>
#include <string>
#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>
#include <vector>
#include "matplotlibcpp.h"

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

namespace plt = matplotlibcpp;
namespace flapping_model{
	const double p = 1.1839;
	const double wing_area = 0.046;
	const double tail_area = 0.018;
	const double C_l0 = 1.332;
	const double D_l0 = 1.713;
	const double D_l1 = -1.639;
	const Vector3d wing_pos{-0.03238, 0.0, 0.0};
	const Vector3d tail_pos{-0.10905, 0.0, 0.02466};
	const double v0 = 2.0;
    const double target_height = 2.0;
	const double f = 10.0; // 扑翼频率
};

void CalF(Model& model, const VectorNd& Q, const VectorNd& QDot, const char* body_name,
			std::vector<SpatialVector>& fext,  
			std::vector<double>& f_x, std::vector<double>& f_z, std::vector<double>& all_angle,
			std::vector<double>& v_x, std::vector<double>& v_z){
	int body_id = model.GetBodyId(body_name);
	double area = flapping_model::wing_area;
	Vector3d body_pos = flapping_model::wing_pos;
	if(body_id == 7){ // body_id of tail is 7
		area = flapping_model::tail_area;
		body_pos = flapping_model::tail_pos;
	}
	// 惯性系下的速度
	Vector3d V_inertial = CalcPointVelocity(model, Q, QDot, body_id, body_pos, true);
	// std::cout << "V_inertial: " << V_inertial.transpose() << std::endl;
	Matrix3d MatWorld2Body = CalcBodyWorldOrientation(model, Q, body_id, true);
	Eigen::Quaterniond q(MatWorld2Body);
	// 从惯性系到body(local)系下的旋转矩阵
	MatWorld2Body = q.normalized().toRotationMatrix();
	// body(local)系下的速度
	Vector3d V_local = MatWorld2Body * V_inertial;// CalcBodyWorldOrientation
	// std::cout << "V_local: " << V_local.transpose() << std::endl;

	Vector3d V_local_proj{V_local[0], 0, V_local[2]};
	double x = V_local[2] / V_local[0];
	double angle_of_attack;
	if(V_local[0] == 0){
		angle_of_attack = atan(INFINITY);
	}
	else{
		angle_of_attack = abs(atan(x));
	}
	Matrix3d rotate_y;
	rotate_y << 0, 0, -1, 0, 1, 0, 1, 0, 0;
	if(V_local[2] > 0.0){
		rotate_y(0, 2) = 1;
		rotate_y(2, 0) = -1;
	}
	// std::cout << rotate_y << std::endl;
	// std::cout << "angle of attack: " << angle_of_attack * 180 / 3.1415 << std::endl;
	double C_l = flapping_model::C_l0 * std::sin(2 * angle_of_attack);
	double C_d = flapping_model::D_l0 - flapping_model::D_l1 * std::cos(2 * angle_of_attack);
	double F_l = 2.0 / 3.0 * C_l * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	double F_d = 2.0 / 3.0 * C_d * flapping_model::p * area * (V_local_proj.norm() * V_local_proj.norm());
	// std::cout << "V: " << V_local_proj.norm() << std::endl;
	// std::cout << "C_l: " << C_l << std::endl;
	// std::cout << "C_d: " << C_d << std::endl;
	// std::cout << "F_l: " << F_l << std::endl;
	// std::cout << "F_d: " << F_d << std::endl;

	// body上的点在惯性系下的位置
	Vector3d pos = CalcBodyToBaseCoordinates(model, Q, body_id, body_pos, true);
	// std::cout << "D: " << D << std::endl;
	Vector3d F = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized());
	// Vector3d F = F_l * (rotate_y * V_local_proj.normalized()) - F_d * (V_local_proj.normalized());
	if(body_id == model.GetBodyId("left_wing")){
		all_angle.push_back(angle_of_attack * 180 / 3.1415);
		Vector3d F_body = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized());
		f_x.push_back(F_body[0]);
		f_z.push_back(F_body[2]);
		v_x.push_back(V_local[0]);
		v_z.push_back(V_local[2]);
	}
	// std::cout << "F1: " << (F_l * (rotate_y * V_local_proj.normalized())).transpose() << std::endl;
	// std::cout << "F2: " << (F_d * (V_local_proj).normalized()).transpose() << std::endl;
	// std::cout << "F1 - F2: " << F_l * (rotate_y * V_local_proj.normalized()) - F_d * (V_local_proj.normalized()) << std::endl;
	// std::cout << "F1 body: " << (F_l *  MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())).transpose() << std::endl;
	// std::cout << "F2 body: " << (F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized())).transpose() << std::endl;
	// std::cout << "F: " << F.transpose() << std::endl;
	Vector3d T = pos.cross(F);
	fext[body_id][0] = T[0];
	fext[body_id][1] = T[1];
	fext[body_id][2] = T[2];
	fext[body_id][3] = F[0];
	fext[body_id][4] = F[1];
	fext[body_id][5] = F[2];
	// std::cout << "fext: " << fext[body_id].transpose() << std::endl;
}

int main (int argc, char* argv[]) {
	rbdl_check_api_version (RBDL_API_VERSION);

	Model* model = new Model();

	if (argc != 2) {
		std::cerr << "Error: not right number of arguments." << std::endl;
		std::cerr << "usage: " << argv[0] << " <model.urdf>" << std::endl;
		exit(-1);
	}

	if (!Addons::URDFReadFromFile (argv[1], model, true, true)) {
		std::cerr << "Error loading model " << argv[1] << std::endl;
		abort();
	}

	std::cout << "Degree of freedom overview:" << std::endl;
	std::cout << Utils::GetModelDOFOverview(*model);

	std::cout << "Model Hierarchy:" << std::endl;
	std::cout << Utils::GetModelHierarchy(*model);
	std::vector<double> wing_f_x, wing_f_z, a_x, alpha_y, tail_alpha_y, a_z, all_theta, body_f_x, body_f_z, all_angle, v_x, v_z;
	for(double theta = 0.0; theta <= 3.14 * 2.0; theta += 3.14 * 2.0 / 50.0){
		// double theta = 0.0;
		all_theta.push_back(theta);
		double theta1 = 0.806 * std::sin(theta) + 0.241;
		double theta3 = -0.806 * std::sin(theta) - 0.241;
		double theta2 = 0.35 * std::sin(theta + 3.1415 / 2.0);
		double theta4 = theta2;
		double theta1_dot = flapping_model::f * 2.0 * 3.1415 * 0.806 * std::cos(theta);
		double theta3_dot = -0.806 * flapping_model::f * 2.0 * 3.1415 * std::cos(theta);
		double theta2_dot = flapping_model::f * 2.0 * 3.1415 * 0.35 * std::cos(theta + 3.1415 / 2.0);
		double theta4_dot = theta2_dot;

		VectorNd Q = VectorNd::Zero (model->q_size);
		Q << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				theta1, theta2, theta3, theta4, -0.0, 0.0;

		VectorNd QDot = VectorNd::Zero (model->qdot_size);
		QDot << 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				theta1_dot, theta2_dot, theta3_dot, theta4_dot, 0.0;

		// // Q[10] = -0.17453;
		// Q << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.35, 0.0, -0.35, 0.0, 0.0;
		// QDot << 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25.0, 0.0, 25.0, 0.0, 0.0;
		// std::cout << "Q: " << Q.transpose() << std::endl;
		// std::cout << "QDot: " << QDot.transpose() << std::endl;
		VectorNd Tau = VectorNd::Zero (model->qdot_size);
		VectorNd QDDot = VectorNd::Zero (model->qdot_size);
		std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
		fext.resize(model->mBodies.size());
		for (unsigned int i = 0; i < fext.size(); ++i) {
			fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
		}
		CalF(*model, Q, QDot, "left_wing", fext, body_f_x, body_f_z, all_angle, v_x, v_z);
		CalF(*model, Q, QDot, "right_wing", fext, body_f_x, body_f_z, all_angle, v_x, v_z);
		CalF(*model, Q, QDot, "tail", fext, body_f_x, body_f_z, all_angle, v_x, v_z);
		wing_f_x.push_back(fext[model->GetBodyId("left_wing")][3]);
		wing_f_z.push_back(fext[model->GetBodyId("left_wing")][5]);
		tail_alpha_y.push_back(fext[model->GetBodyId("tail")][1]);

		// fext[model->GetBodyId("left_wing")][5] = 3.2;
		
		// fext[model->GetBodyId("body")][5] = 3.2;
		// fext[model->GetBodyId("right_wing")][5] = 3.2;	
		ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
		// ForwardDynamics (*model, Q, QDot, Tau, QDDot); 
		// std::cout << "after QDDot: " << QDDot.transpose() << std::endl;

		// InverseDynamics (*model, Q, QDot, QDDot, Tau);
		a_x.push_back(QDDot[0]);
		a_z.push_back(QDDot[2]);
		alpha_y.push_back(QDDot[4]);
		// std::cout << "QDDOT: " << QDDot.transpose() << std::endl;
	}
	

	delete model;
	plt::subplot(2,3,1);
	plt::plot(all_theta, wing_f_x, {{"label","wing_f_x"}});
	plt::plot(all_theta, wing_f_z, {{"label","wing_f_z"}});
	plt::plot(all_theta, body_f_x, {{"label","body_f_x"}});
	plt::plot(all_theta, body_f_z, {{"label","body_f_z"}});
	plt::legend();
	plt::subplot(2,3,2);
	plt::plot(all_theta, a_x, {{"label","a_x"}});
	plt::plot(all_theta, a_z, {{"label","a_z"}});
	plt::legend();
	plt::subplot(2,3,3);
	plt::plot(all_theta, all_angle, {{"label","angle_of_attack"}});
	plt::legend();
	plt::subplot(2,3,4);
	plt::plot(all_theta, v_x, {{"label","v_local_x"}});
	plt::plot(all_theta, v_z, {{"label","v_local_z"}});
	plt::legend();
	plt::subplot(2,3,5);
	plt::plot(all_theta, alpha_y, {{"label","alpha_y"}});
	plt::plot(all_theta, tail_alpha_y, {{"label","tail_alpha_y"}});
	plt::legend();
    plt::show();

	return 0;
}

