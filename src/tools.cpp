#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// This code is more or less the same as explained on the conferences.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() == 0 || ground_truth.size() == 0){
      cout << "ERROR!! Vector size is zero." << endl;
      return rmse;
    }
	if(estimations.size() != ground_truth.size()){
      cout << "ERROR!! Vector size is not equal." << endl;
      return rmse;
    }

    for(int i=0; i < estimations.size(); ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      rmse += residual;
    }

	//Calculate mean and square root of rmse and return it
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

// This code is more or less the same as explained on the conferences.
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
    double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double denom1 = px*px+py*py;
	double term2 = sqrt(denom1);
	double term3 = (denom1*term2);

	//check division by zero
	if(denom1 == 0.0){
		cout << "ERROR!! Division by Zero." << endl;
		return Hj;
	}

	//compute Jacobian matrix
	Hj << (px/term2), (py/term2), 0, 0,
		  -(py/denom1), (px/denom1), 0, 0,
		  py*(vx*py - vy*px)/term3, px*(px*vy - py*vx)/term3, px/term2, py/term2;

	return Hj;
}    




	
	
	