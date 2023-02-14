/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <cmath>
#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

namespace flapping_model{
	const double p = 1.1839;
	const double area = 0.1;
	const double C_l0 = 1.0;
	const double D_l0 = 1.0;
	const double D_l1 = 1.0;
};


using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

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

	VectorNd Q = VectorNd::Zero (model->q_size);
	VectorNd QDot = VectorNd::Zero (model->qdot_size);
	// QDot << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
	VectorNd Tau = VectorNd::Zero (model->qdot_size);
	VectorNd QDDot = VectorNd::Zero (model->qdot_size);
	Matrix3d rotate_y;
	rotate_y << 0, 0, -1, 0, 1, 0, 1, 0, 0;
	auto V_l_inertial = CalcPointVelocity(*model, Q, QDot, model->GetBodyId("link1"), {1, 1, 0});
	auto V_r_inertial = CalcPointVelocity(*model, Q, QDot, model->GetBodyId("link2"), {0, 0, 0}, false);
	auto V_l_local = CalcBodyWorldOrientation(*model, Q, model->GetBodyId("link1"), false) * V_l_inertial;
	auto V_r_local = CalcBodyWorldOrientation(*model, Q, model->GetBodyId("link2"), false) * V_r_inertial;
	std::cout << "V_l_inertial : " << V_l_inertial << std::endl;
	std::cout << "V_l_local : " << V_l_local << std::endl;
	V_l_inertial.normalize();
	double angle_of_attack;
	if(V_l_local[0] == 0){
		double angle_of_attack = atan(INFINITY);
	}
	else{
		double angle_of_attack = atan(V_l_local[2] / V_l_local[0]);
	}
	
	double C_l = flapping_model::C_l0 * sin(2 * angle_of_attack);
	double C_d = flapping_model::D_l0 - flapping_model::D_l1 * cos(2 * angle_of_attack);


	double F_l_l = 2.0 / 3.0 * C_l * flapping_model::p * flapping_model::area * (V_l_local[0] * V_l_local[0] + V_l_local[2] * V_l_local[2]);
	double F_l_d = 2.0 / 3.0 * C_d * flapping_model::p * flapping_model::area * (V_l_local[0] * V_l_local[0] + V_l_local[2] * V_l_local[2]);
	std::cout << "F_l_l : " << F_l_l << std::endl;
	auto F_l = F_l_l * rotate_y * V_l_inertial - F_l_d * V_l_inertial;
	std::vector<RigidBodyDynamics::Math::SpatialVector> fext;
	fext.resize(model->mBodies.size());
	for (unsigned int i = 0; i < fext.size(); ++i) {
		fext[i]=RigidBodyDynamics::Math::SpatialVector::Zero();
	}
	std::cout << "F_l : " << F_l << std::endl;
	// fext[3][3] = F_l[0];
	// fext[3][4] = F_l[1];
	// fext[3][5] = F_l[2];
	std::cout << "Forward Dynamics with q, qdot, tau set to zero:" << std::endl;
	ForwardDynamics (*model, Q, QDot, Tau, QDDot, &fext);
	// InverseDynamics (*model, Q, QDot, QDDot, Tau);
	std::cout << QDDot.transpose() << std::endl;

	delete model;

	return 0;
}

