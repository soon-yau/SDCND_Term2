#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;
Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	//accumulate squared residuals
	  VectorXd rmse(4);
    rmse << 0,0,0,0;  

    for(unsigned int i=0; i < estimations.size(); ++i){

      VectorXd residual = estimations[i] - ground_truth[i];

      //coefficient-wise multiplication
      residual = residual.array()*residual.array();
      rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj_ = MatrixXd(3, 4);

  float px=x_state(0);
  float py=x_state(1);
  float vx=x_state(2);
  float vy=x_state(3);
  float px2py2=px*px+py*py;
  float sqrt_px2py2=sqrt(px2py2);
  float px2py2_3_2=pow(px2py2,1.5);
  
  // if divide by zero, return empty Jacobian matrix
  if (px2py2<0.0001) return Hj_;

  float px_sqrt_px2py2=px/sqrt_px2py2;
  float py_sqrt_px2py2=py/sqrt_px2py2;
  
  Hj_(0,0)=px_sqrt_px2py2;
  Hj_(0,1)=py_sqrt_px2py2;
  Hj_(0,2)=0;
  Hj_(0,3)=0;

  Hj_(1,0)=-py/px2py2;
  Hj_(1,1)=px/px2py2;
  Hj_(1,2)=0;
  Hj_(1,3)=0;

  Hj_(2,0)=py*(vx*py-vy*px)/px2py2_3_2;
  Hj_(2,1)=px*(vy*px-vx*py)/px2py2_3_2;
  Hj_(2,2)=px_sqrt_px2py2;
  Hj_(2,3)=py_sqrt_px2py2;

  return Hj_;
}
