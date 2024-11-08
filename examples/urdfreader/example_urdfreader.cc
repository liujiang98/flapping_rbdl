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
	const double tail_area = 0.018; //0.03845369   0.018
	const double C_l0 = 1.332;
	const double D_l0 = 1.713;
	const double D_l1 = -1.639;
	const Vector3d wing_pos{-0.03238, 0.0, 0.0};
	const Vector3d tail_pos{-0.10905, 0.0, 0.02466}; //-0.2635, 0, 0.0268   -0.10905, 0.0, 0.02466
	const double v0 = 2.0;
    const double target_height = 2.0;
	const double f = 12.0; // 扑翼频率
};

void CalF(Model& model, const VectorNd& Q, const VectorNd& QDot, const char* body_name,
			std::vector<SpatialVector>& fext,  
			std::vector<double>& f_x, std::vector<double>& f_z, std::vector<double>& all_angle,
			std::vector<double>& v_x, std::vector<double>& v_z, std::vector<double>& cd_wing,
			std::vector<double>& cl_wing, std::vector<double>& cd_tail, std::vector<double>& cl_tail,
			std::vector<double>& wing_l_t_y, std::vector<double>& wing_d_t_y, std::vector<double>& wing_pos_z){
	int body_id = model.GetBodyId(body_name);
	double area = flapping_model::wing_area;
	Vector3d body_pos = flapping_model::wing_pos;
	if(body_id == 7){ // body_id of tail is 7
		area = flapping_model::tail_area;
		body_pos = flapping_model::tail_pos;
	}
	// 惯性系下的速度
	Vector3d V_inertial = CalcPointVelocity(model, Q, QDot, body_id, body_pos, true);
	std::cout << "V_inertial: " << V_inertial.transpose() << std::endl;
	Matrix3d MatWorld2Body = CalcBodyWorldOrientation(model, Q, body_id, true);
	Eigen::Quaterniond q(MatWorld2Body);
	// 从惯性系到body(local)系下的旋转矩阵
	MatWorld2Body = q.normalized().toRotationMatrix();
	// body(local)系下的速度
	Vector3d V_local = MatWorld2Body * V_inertial;// CalcBodyWorldOrientation
	std::cout << "V_local: " << V_local.transpose() << std::endl;

	Vector3d V_local_proj{V_local[0], 0, V_local[2]};
	double x;
	double angle_of_attack;
	if(V_local[0] == 0){
		angle_of_attack = atan(INFINITY);
	}
	else{
		x = V_local[2] / V_local[0];
	}
	angle_of_attack = atan(x);
	Matrix3d rotate_y;
	rotate_y << 0, 0, 1, 0, 1, 0, -1, 0, 0;
	// if(V_local[2] > 0.0){
	// 	rotate_y(0, 2) = 1;
	// 	rotate_y(2, 0) = -1;
	// }
	
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
	// if(body_id == 7){ // body_id of tail is 7
	// 	// pos[0] = 5.5 * pos[0];
	// 	F_l *= 12.0;
	// 	// F_d *= 5.0;
	// }
	
	std::cout << "pos: " << pos.transpose() << std::endl;
	// 确定力在惯性系下的方向
	Vector3d F = F_l * MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())
			- F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized());
	// Vector3d F1 = F_l * (rotate_y * V_local_proj.normalized()) - F_d * (V_local_proj.normalized());
	if(body_id == model.GetBodyId("left_wing")){
		all_angle.push_back(angle_of_attack * 180 / 3.1415);
		f_x.push_back(F[0]);
		f_z.push_back(F[2]);
		v_x.push_back(V_local[0]);
		v_z.push_back(V_local[2]);
		cd_wing.push_back(C_d);
		cl_wing.push_back(C_l);
		wing_l_t_y.push_back(-pos[0] * F[2]);
		wing_d_t_y.push_back(pos[2] * F[0]);
		wing_pos_z.push_back(pos[2]);
	}
	// if(body_id == model.GetBodyId("tail")){
	// 	cd_tail.push_back(C_d);
	// 	cl_tail.push_back(C_l);
	// }
	// std::cout << "F1: " << (F_l * (rotate_y * V_local_proj.normalized())).transpose() << std::endl;
	// std::cout << "F2: " << (F_d * (V_local_proj).normalized()).transpose() << std::endl;
	// std::cout << "F1 - F2: " << F_l * (rotate_y * V_local_proj.normalized()) - F_d * (V_local_proj.normalized()) << std::endl;
	// std::cout << "F1 body: " << (F_l *  MatWorld2Body.transpose() * (rotate_y * V_local_proj.normalized())).transpose() << std::endl;
	// std::cout << "F2 body: " << (F_d * ((MatWorld2Body.transpose() * V_local_proj).normalized())).transpose() << std::endl;
	std::cout << "F: " << F.transpose() << std::endl;
	Vector3d T = pos.cross(F);
	std::cout << "T: " << T.transpose() << std::endl;
	// std::cout << "D * F: " << (D * F).transpose() << std::endl;
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

	if (!Addons::URDFReadFromFile (argv[1], model, true, false)) {
		std::cerr << "Error loading model " << argv[1] << std::endl;
		abort();
	}

	// std::cout << "Degree of freedom overview:" << std::endl;
	// std::cout << Utils::GetModelDOFOverview(*model);

	// std::cout << "Model Hierarchy:" << std::endl;
	// std::cout << Utils::GetModelHierarchy(*model);
	
	std::vector<double> cd_wing, cl_wing, cd_tail, cl_tail, wing_l_t_y, wing_d_t_y, 
						wing_f_x, wing_f_z, a_x, alpha_y, wing_t_y, tail_t_y, a_z, 
						all_theta, body_f_x, body_f_z, all_angle, v_x, v_z, wing_pos_z, tail_theta, tail_t;
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
				theta1, theta2, theta3, theta4, -0.58, 1;
		// Q << 0.0, 0.0, 0.0, 0.0, -0.087, 0.0, 
		// 		theta1, theta2, theta3, theta4, 1.0, 0.9962;

		VectorNd QDot = VectorNd::Zero (model->qdot_size);
		QDot << 3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
		CalF(*model, Q, QDot, "left_wing", fext, body_f_x, body_f_z, all_angle, v_x, 
				v_z, cd_wing, cl_wing, cd_tail, cl_tail, wing_l_t_y, wing_d_t_y, wing_pos_z);
		CalF(*model, Q, QDot, "right_wing", fext, body_f_x, body_f_z, all_angle, v_x, 
				v_z, cd_wing, cl_wing, cd_tail, cl_tail, wing_l_t_y, wing_d_t_y, wing_pos_z);
		CalF(*model, Q, QDot, "tail", fext, body_f_x, body_f_z, all_angle, v_x, 
				v_z, cd_wing, cl_wing, cd_tail, cl_tail, wing_l_t_y, wing_d_t_y, wing_pos_z);
		std::cout << "///////////////////////////////" << std::endl;
		// fext[model->GetBodyId("tail")][1] = 0.0;
		// fext[model->GetBodyId("tail")][3] = 0.0;
		// fext[model->GetBodyId("tail")][5] = 0.0;
		// fext[model->GetBodyId("left_wing")][1] = 0.0;
		// fext[model->GetBodyId("left_wing")][3] = 0.0;
		// fext[model->GetBodyId("left_wing")][5] = 0.0;
		// fext[model->GetBodyId("right_wing")][1] = 0.0;
		// fext[model->GetBodyId("right_wing")][3] = 0.0;
		// fext[model->GetBodyId("right_wing")][5] = 0.0;
		wing_f_x.push_back(fext[model->GetBodyId("left_wing")][3]);
		wing_f_z.push_back(fext[model->GetBodyId("left_wing")][5]);
		tail_t_y.push_back(fext[model->GetBodyId("tail")][1]);
		wing_t_y.push_back(fext[model->GetBodyId("left_wing")][1]);
		

		// fext[model->GetBodyId("left_wing")][1] = 3.2;
		
		// fext[model->GetBodyId("body")][5] = 3.2;
		// fext[model->GetBodyId("right_wing")][1] = 3.2;
		ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
		// for(int i = 1; i < model->mBodies.size(); i++){
		// 	std::cout << "External force (" << i << ") = " << model->X_base[i].toMatrixAdjoint() * (fext)[i] << std::endl;
		// 	Vector3d T{fext[i][0], fext[i][1], fext[i][2]};
		// 	Vector3d F{fext[i][3], fext[i][4], fext[i][5]};
		// 	// std::cout << "for T: " << T.transpose() << std::endl;
		// 	std::cout << "T 111: " << model->X_base[i].E * T << std::endl;
		// 	std::cout << "X base E: " << model->X_base[i].E << std::endl;
		// 	std::cout << "X base r: " << model->X_base[i].r << std::endl;
		// 	Matrix3d _Erx =
		// 		model->X_base[i].E * Matrix3d (
		// 		0., -model->X_base[i].r[2], model->X_base[i].r[1],
		// 		model->X_base[i].r[2], 0., -model->X_base[i].r[0],
		// 		-model->X_base[i].r[1], model->X_base[i].r[0], 0.
		// 		);
		// 	std::cout << "T 222: " << _Erx * F << std::endl;
		// }
		// ForwardDynamics (*model, Q, QDot, Tau, QDDot); 
		// std::cout << "after QDDot: " << QDDot.transpose() << std::endl;

		// InverseDynamics (*model, Q, QDot, QDDot, Tau);
		a_x.push_back(QDDot[0]);
		a_z.push_back(QDDot[2]);
		alpha_y.push_back(QDDot[4]);
		// std::cout << "QDDOT: " << QDDot.transpose() << std::endl;
	}

	for(double theta = -1.5; theta <= 1.5; theta += 0.05){
		tail_theta.push_back(theta);
		VectorNd Q = VectorNd::Zero (model->q_size);
		Q << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.0, 0.0, 0.0, 0.0, theta, 1;
		// Q << 0.0, 0.0, 0.0, 0.0, -0.087, 0.0, 
		// 		theta1, theta2, theta3, theta4, 1.0, 0.9962;

		VectorNd QDot = VectorNd::Zero (model->qdot_size);
		QDot << 3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 0.0, 0.0, 0.0;

		VectorNd Tau = VectorNd::Zero (model->qdot_size);
		VectorNd QDDot = VectorNd::Zero (model->qdot_size);
		std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
		fext.resize(model->mBodies.size());
		for (unsigned int i = 0; i < fext.size(); ++i) {
			fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
		}
		CalF(*model, Q, QDot, "tail", fext, body_f_x, body_f_z, all_angle, v_x, 
				v_z, cd_wing, cl_wing, cd_tail, cl_tail, wing_l_t_y, wing_d_t_y, wing_pos_z);
		tail_t.push_back(fext[model->GetBodyId("tail")][1]);
		// ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
		
	}
	

	delete model;
	plt::subplot(2,4,1);
	plt::plot(all_theta, wing_f_x, {{"label","wing_f_x"}});
	plt::plot(all_theta, wing_f_z, {{"label","wing_f_z"}});
	plt::plot(all_theta, body_f_x, {{"label","body_f_x"}});
	plt::plot(all_theta, body_f_z, {{"label","body_f_z"}});
	plt::legend();
	plt::subplot(2,4,2);
	plt::plot(all_theta, a_x, {{"label","a_x"}});
	plt::plot(all_theta, a_z, {{"label","a_z"}});
	plt::legend();
	plt::subplot(2,4,3);
	plt::plot(all_theta, all_angle, {{"label","angle_of_attack"}});
	plt::legend();
	plt::subplot(2,4,4);
	plt::plot(all_theta, v_x, {{"label","v_local_x"}});
	plt::plot(all_theta, v_z, {{"label","v_local_z"}});
	plt::legend();
	plt::subplot(2,4,5);
	plt::plot(all_theta, alpha_y, {{"label","alpha_y"}});
	plt::legend();
	plt::subplot(2,4,6);
	plt::plot(all_theta, wing_t_y, {{"label","wing_t_y"}});
	plt::plot(all_theta, tail_t_y, {{"label","tail_t_y"}});
	plt::plot(all_theta, wing_l_t_y, {{"label","wing_l_t_y"}});
	plt::plot(all_theta, wing_d_t_y, {{"label","wing_d_t_y"}});
	plt::legend();
	plt::subplot(2,4,7);
	plt::plot(all_theta, cd_wing, {{"label","cd_wing"}});
	plt::plot(all_theta, cl_wing, {{"label","cl_wing"}});
	// plt::plot(all_theta, cd_tail, {{"label","cd_tail"}});
	// plt::plot(all_theta, cl_tail, {{"label","cl_tail"}});
	plt::legend();
	plt::subplot(2,4,8);
	plt::plot(tail_theta, tail_t, {{"label","tail_t"}});
	plt::legend();
    plt::show();

	return 0;
}

