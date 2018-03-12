#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Random(5, 5);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.52;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  ///* State dimension
  n_x_=5;

  ///* Augmented state dimension
  n_aug_=7;

  ///* Sigma point spreading parameter
  lambda_=3-n_x_;

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0)=lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++)
  {
      weights_(i)=0.5/(lambda_+n_aug_);
  }
  
  // predicted sigma 
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);  

  R_radar = MatrixXd(3,3);
  R_radar << std_radr_*std_radr_, 0, 0,
             0,  std_radphi_*std_radphi_, 0,
             0, 0, std_radrd_*std_radrd_;  
        
  R_laser = MatrixXd(2,2);
  R_laser << std_laspx_*std_laspx_, 0,
             0, std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_)
  {
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {  
      double rho=meas_package.raw_measurements_[0];
      double phi=meas_package.raw_measurements_[1];
      double rho_dot=meas_package.raw_measurements_[2];

      double px = rho*cos(phi);
      double py = rho*sin(phi);

      x_ << px, py, 0, 0, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px=meas_package.raw_measurements_[0];
      double py=meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    }
    is_initialized_=true;

    return;
  }

  double delta_t=(meas_package.timestamp_-time_us_)/1000000.0;
  //cout<<delta_t<<endl;
  time_us_ = meas_package.timestamp_;
  
  Prediction(delta_t);
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {  
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Create augmentation
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_)=x_;
  x_aug(5)=0;
  x_aug(6)=0;
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6)=std_yawdd_*std_yawdd_;

  // Generate Sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd A = P_aug.llt().matrixL();
  
  Xsig_aug.col(0)=x_aug;  
  for (int i=1; i<=n_aug_; i++)
  {
      Xsig_aug.col(i)=x_aug+sqrt(lambda_+n_aug_)*A.col(i-1);
      Xsig_aug.col(i+n_aug_)=x_aug-sqrt(lambda_+n_aug_)*A.col(i-1);
  }
  
  // Predict Sigma points
  VectorXd noise=VectorXd(n_x_);
  VectorXd pred=VectorXd(n_x_);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      double px=Xsig_aug(0, i);
      double py=Xsig_aug(1, i);
      double v=Xsig_aug(2, i);
      double phi=Xsig_aug(3, i);
      double phi_dot=Xsig_aug(4, i);
      double noise_a=Xsig_aug(5, i);
      double noise_phi_ddot=Xsig_aug(6, i);

      noise(0)=0.5*delta_t*delta_t*cos(phi)*noise_a;
      noise(1)=0.5*delta_t*delta_t*sin(phi)*noise_a;
      noise(2)=delta_t*noise_a;
      noise(3)=0.5*delta_t*delta_t*noise_phi_ddot;
      noise(4)=delta_t*noise_phi_ddot;
      
      if (fabs(phi_dot)<0.001)
      {
          pred(0)=v*cos(phi)*delta_t;
          pred(1)=v*sin(phi)*delta_t;
          pred(2)=0;
          pred(3)=phi_dot*delta_t;
          pred(4)=0;          
      }
      
      else
      {
          
          pred(0)=(v/phi_dot)*(sin(phi+phi_dot*delta_t)-sin(phi));
          pred(1)=(v/phi_dot)*(-cos(phi+phi_dot*delta_t)+cos(phi));
          pred(2)=0;
          pred(3)=phi_dot*delta_t;
          pred(4)=0;
      }
      //VectorXd X= VectorXd(n_x_);
      //X<< px, py, v, phi, phi_dot;
      //Xsig_pred_.col(i)=X+pred+noise;      
      Xsig_pred_.col(i)=Xsig_aug.block<5,1>(0,i)+pred+noise;      
  }  

  // Calculate mean and covariance
  //predict state mean
  x_.fill(0.0);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      x_=x_+weights_(i)*Xsig_pred_.col(i);
  }
  
  //predict state covariance matrix
  P_.fill(0.0);
  VectorXd vX_x = VectorXd(n_x_);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      vX_x=Xsig_pred_.col(i)-x_;
      
      while (vX_x(3)> M_PI) vX_x(3)-=2.*M_PI;
      while (vX_x(3)<-M_PI) vX_x(3)+=2.*M_PI;     
      P_=P_+weights_(i)*vX_x*vX_x.transpose();
  }
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //create  vector for incoming radar measurement
  int n_z = 2;
  VectorXd z = VectorXd(n_z);
  z <<  meas_package.raw_measurements_[0],  // px
        meas_package.raw_measurements_[1];  // py

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);  
  z_pred.fill(0.0);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      // copy from Xsig_pred_ (px and py)
      Zsig.col(i)=Xsig_pred_.block<2,1>(0,i);

      z_pred=z_pred+weights_(i)*Zsig.col(i);
  }

  //calculate measurement covariance matrix S      
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  VectorXd vZ_z = VectorXd(n_z);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {   
      //residual
      vZ_z=Zsig.col(i)-z_pred;
      S=S+weights_(i)*vZ_z*vZ_z.transpose();
  }
    
  S=S+R_laser;

  //calculate cross correlation matrix
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd Z_z = VectorXd(n_z);
  VectorXd X_x = VectorXd(n_x_); 
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++)
  {
      X_x = Xsig_pred_.col(i)-x_;
      Z_z = Zsig.col(i)-z_pred;
      Tc=Tc+weights_(i)*X_x*Z_z.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_z, n_z);
  K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_laser_=z_diff.transpose()*S.inverse()*z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  double px, py, v, phi;
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      px=Xsig_pred_(0,i);
      py=Xsig_pred_(1,i);
      v=Xsig_pred_(2,i);
      phi=Xsig_pred_(3,i);
      
      Zsig(0,i)=sqrt(px*px+py*py);
      Zsig(1,i)=atan2(py,px);
      if (fabs(Zsig(0,i))<0.0001)
        Zsig(2,i)=(px*cos(phi)*v+py*sin(phi)*v)/0.0001;
      else
        Zsig(2,i)=(px*cos(phi)*v+py*sin(phi)*v)/Zsig(0,i);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      z_pred=z_pred+weights_(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S      
  VectorXd vZ_z = VectorXd(n_z);

  S.fill(0.0);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {   
      //residual
      vZ_z=Zsig.col(i)-z_pred;
      //angle normalization
      while (vZ_z(1)> M_PI) vZ_z(1)-=2.*M_PI;
      while (vZ_z(1)<-M_PI) vZ_z(1)+=2.*M_PI;
      S=S+weights_(i)*vZ_z*vZ_z.transpose();
  }
    
  S=S+R_radar;

  //create  vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<  meas_package.raw_measurements_[0],  // rho
        meas_package.raw_measurements_[1],  // phi
        meas_package.raw_measurements_[2];  // rho_dot

  //calculate cross correlation matrix
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd Z_z = VectorXd(n_z);
  VectorXd X_x = VectorXd(n_x_); 
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++)
  {
      X_x = Xsig_pred_.col(i)-x_;
      while (X_x(3)> M_PI) X_x(3)-=2.*M_PI;
      while (X_x(3)<-M_PI) X_x(3)+=2.*M_PI;           
      Z_z = Zsig.col(i)-z_pred;
      while (Z_z(1)> M_PI) Z_z(1)-=2.*M_PI;
      while (Z_z(1)<-M_PI) Z_z(1)+=2.*M_PI;     
      Tc=Tc+weights_(i)*X_x*Z_z.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_radar_=z_diff.transpose()*S.inverse()*z_diff;
}
