#include <iostream>
#include "tools.h"


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	// Calculate the RMSE
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size() + 1
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	vector<VectorXd> estm;
	for (int i=0; i < estimations.size(); ++i) {
    	VectorXd converted(4);
    	converted << estimations[i][0],                           // px
                     estimations[i][1],                           // py
                     cos(estimations[i][3])*estimations[i][2],    // vx
                     sin(estimations[i][3])*estimations[i][2];    // vy
    	estm.push_back(converted);
  	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estm.size(); ++i){

		VectorXd residual = estm[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estm.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;

}
