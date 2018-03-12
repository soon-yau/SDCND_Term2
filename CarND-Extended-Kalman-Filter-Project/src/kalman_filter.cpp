#include <iostream>
#include <cmath>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  cout<<"Init"<<endl;  
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  //cout<<"Predict"<<endl;
  x_ = F_ * x_;  
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  //cout<<"Update"<<endl;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //cout<<"Update"<<endl;
  float px=x_(0);
  float py=x_(1);
  float vx=x_(2);
  float vy=x_(3);
  
  VectorXd z_pred=VectorXd(3); 
  z_pred(0)=sqrt(px*px+py*py);
  z_pred(1)=atan2(py,px);
  float denominator=z_pred(0)<0.001?0.001:z_pred(0);
  z_pred(2)=(px*vx+py*vy)/denominator;

  VectorXd y = z - z_pred;
  Tools tools;
  MatrixXd Hj_= tools.CalculateJacobian(x_);
  MatrixXd Ht = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  //cout<<"\nK:"<<K<<endl;
  //cout<<"\nY:"<<y<<endl;
  x_ = x_ + (K * y);
  //cout<<"\nx_:"<<x_<<endl;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;    

}
