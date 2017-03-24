#include <iostream>
#include "ukf.h"

using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // previous timestamp
  previous_timestamp_ = 0;

  // initial state vector
  //set the state vector with CTRV model: px, py, speed, yaw angle, yaw angle rate
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ <<    0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  // State dimension
  n_x_ = 5;

  //Augmented state dimension
  n_aug_ = 7;

  //set radar measurement dimension, radar can measure r, phi, and r_dot
  n_radar_ = 3;

  //set laser measurement dimension, laser can measure px and py
  n_laser_ = 2;

  //define augment spreading parameter
  lambda_ = 3 - n_aug_;

  //create augment sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  //matrix for sigma points in radar measurement space
  radar_sig_ = MatrixXd(n_radar_, 2 * n_aug_ + 1);

  //matrix for sigma points in laser measurement space
  laser_sig_ = MatrixXd(n_laser_, 2 * n_aug_ + 1);

  //mean predicted radar measurement
  radar_pred_ = VectorXd(n_radar_);

  //mean predicted laser measurement
  laser_pred_ = VectorXd(n_laser_);

  //radar measurement covariance matrix
  radar_S_ = MatrixXd(n_radar_,n_radar_);

  //laser measurement covariance matrix
  laser_S_ = MatrixXd(n_laser_,n_laser_);

  // Radar NIS
  radar_nis_ = 0;

  // laser NIS
  laser_nis_ = 0;
}

UKF::~UKF() {}

/*
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  //calculate square root of P_
  MatrixXd A_ = P_.llt().matrixL();

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_ + n_x_) * A_.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A_.col(i);
  }

  //write result
  *Xsig_out = Xsig;

}
*/

void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

}

void UKF::SigmaPointPrediction(double delta_t) {
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}

void UKF::PredictMeanAndCovariance() {
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  //cout << "here 5" <<endl;
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }
  //cout << "here 6" <<endl;
  //predicted state covariance matrix
  P_.fill(0.0);
  cout << "Xsig_pred_" <<Xsig_pred_ <<endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //cout << "here 7" <<endl;
    //cout << "x_diff(3)" <<x_diff(3) <<endl;
    //cout << "x_" <<x_ <<endl;


    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //cout << "here 8" <<endl;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    //cout << "here 9" <<endl;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
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

  AugmentedSigmaPoints();
  //cout << "here 3" <<endl;
  SigmaPointPrediction(delta_t);
  //cout << "here 4" <<endl;
  PredictMeanAndCovariance();


}

void UKF::PredictRadarMeasurement() {
  //set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // radar measurement model
    radar_sig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    radar_sig_(1,i) = atan2(p_y,p_x);                                 //phi
    radar_sig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //radar mean predicted measurement
  radar_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      radar_pred_ = radar_pred_ + weights_(i) * radar_sig_.col(i);
  }

  //radar measurement covariance matrix
  radar_S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = radar_sig_.col(i) - radar_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    radar_S_ = radar_S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_radar_,n_radar_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  radar_S_ = radar_S_ + R;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_radar_);
  double r = meas_pack.raw_measurements_[0];
  double phi = meas_pack.raw_measurements_[1];
  double rd = meas_pack.raw_measurements_[2];
  z << 
      r,
      phi,
      rd;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_radar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = radar_sig_.col(i) - radar_pred_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * radar_S_.inverse();

  //residual
  VectorXd z_diff = z - radar_pred_;

  //Calculate radar NIS
  radar_nis_ = z_diff.transpose()*radar_S_.inverse()*z_diff;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*radar_S_*K.transpose();    

}

void UKF::PredictLaserMeasurement() {
  //set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // laser measurement model
    laser_sig_(0,i) = p_x; //px
    laser_sig_(1,i) = p_y; //py

  }

  //laser mean predicted measurement
  laser_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      laser_pred_ = laser_pred_ + weights_(i) * laser_sig_.col(i);
  }

  //laser measurement covariance matrix
  laser_S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = laser_sig_.col(i) - laser_pred_;

    laser_S_ = laser_S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_laser_,n_laser_);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  laser_S_ = laser_S_ + R;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_laser_);
  double p_x = meas_pack.raw_measurements_[0];
  double p_y = meas_pack.raw_measurements_[1];
  z << 
      p_x,
      p_y;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_laser_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = laser_sig_.col(i) - laser_pred_;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * laser_S_.inverse();

  //residual
  VectorXd z_diff = z - laser_pred_;

  //Calculate laser NIS
  laser_nis_ = z_diff.transpose()*laser_S_.inverse()*z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*laser_S_*K.transpose();
}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
   ****************************************************************************/

  //cout << "measurement: " << measurement_pack.raw_measurements_ <<endl;
  if (!is_initialized_) {

    // first measurement
    //cout << "UKF: " << endl;

    double px;
    double py;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      px = ro * cos(phi);
      py = ro * sin(phi);
      x_ << px, py, 0, 0, 0;
      //cout << "radar init" <<endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      //set the state with the initial location and zero speed, zero yaw angle, zero yaw angle rate
      x_ << px, py, 0, 0, 0;
      //cout << "laser init" <<endl;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    //cout << "init finished" <<endl;
    //cout << "x_" << x_ <<endl;
    //cout << "P_" << P_ <<endl;
    
    return;
  }

  // Skip 0's measurements
  if (abs(measurement_pack.raw_measurements_[0]) < 0.000001 or abs(measurement_pack.raw_measurements_[1]) < 0.000001) {
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //cout << "here 1" <<endl;
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  //cout << "here 2" <<endl;
  cout << "dt" <<dt <<endl;
  Prediction(dt);



  /*****************************************************************************
   *  Update
      * Use the sensor type to perform the update step.
      * Update the state and covariance matrices.
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if (use_radar_) {
      PredictRadarMeasurement();
      UpdateRadar(measurement_pack);
      //cout << "radar update" <<endl;
    }

  } else {
    // Laser updates
    if (use_laser_) {
      PredictLaserMeasurement();
      UpdateLidar(measurement_pack);
      //cout << "laser update" <<endl;
    }
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}
