#include <UnitTest++.h>

#include <iostream>

#include "mathutils.h"
#include "Logging.h"

#include "Model.h"
#include "Contacts.h"
#include "Dynamics.h"
#include "Kinematics.h"

using namespace std;
using namespace SpatialAlgebra;

const double TEST_PREC = 1.0e-14;

struct ContactsFixture {
	ContactsFixture () {
		ClearLogOutput();
		model = new Model;
		model->Init();

		model->gravity.set (0., -9.81, 0.);

		/* Basically a model like this, where X are the Center of Masses
		 * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
		 *
		 *                X
		 *                *
		 *              _/
		 *            _/  (-Z)
		 *      Z    /
		 *      *---* 
		 *      |
		 *      |
		 *  Z   |
		 *  O---*
		 *      Y
		 */

		body_a = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
		joint_a = Joint(
				JointTypeRevolute,
				Vector3d (0., 0., 1.)
				);

		body_a_id = model->AddBody(0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);

		body_b = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
		joint_b = Joint (
				JointTypeRevolute,
				Vector3d (0., 1., 0.)
				);

		body_b_id = model->AddBody(body_a_id, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

		body_c = Body (1., Vector3d (0., 0., 1.), Vector3d (1., 1., 1.));
		joint_c = Joint (
				JointTypeRevolute,
				Vector3d (0., 0., 1.)
				);

		body_c_id = model->AddBody(body_b_id, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

		body_d = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
		joint_c = Joint (
				JointTypeRevolute,
				Vector3d (1., 0., 0.)
				);

		body_d_id = model->AddBody(body_c_id, Xtrans(Vector3d(0., 0., -1.)), joint_c, body_d);

		Q = cmlVector(4);
		QDot = cmlVector(4);
		QDDot = cmlVector(4);
		Tau = cmlVector(4);

		Q.zero();
		QDot.zero();
		QDDot.zero();
		Tau.zero();

		body_id = body_c_id;
		contact_point.set (1., 0., 0.);
		contact_normal.set (0., 1., 0.);

		ClearLogOutput();
	}
	
	~ContactsFixture () {
		delete model;
	}
	Model *model;

	unsigned int body_a_id, body_b_id, body_c_id, body_d_id;
	Body body_a, body_b, body_c, body_d;
	Joint joint_a, joint_b, joint_c, joint_d;

	cmlVector Q;
	cmlVector QDot;
	cmlVector QDDot;
	cmlVector Tau;

	unsigned int body_id;
	Vector3d contact_point;
	Vector3d contact_normal;
};

TEST_FIXTURE(ContactsFixture, TestContactSimple) {
	model->AddContact(body_id, contact_point, contact_normal);

	Vector3d point_accel;
	{
//		SUPPRESS_LOGGING;
			_NoLogging nolog;
		CalcPointAcceleration (*model, Q, QDot, QDDot, body_id, contact_point, point_accel);
	}
	cout << "point accel pre  = " << point_accel << endl;

	ForwardDynamicsContacts(*model, Q, QDot, Tau, QDDot);
	{
		SUPPRESS_LOGGING
		_NoLogging nolog;
		CalcPointAcceleration (*model, Q, QDot, QDDot, body_id, contact_point, point_accel);
	}
	cout << LogOutput.str() << endl;
	cout << "point accel post = " << point_accel << endl;
}