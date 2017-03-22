#include <iostream>
#include "ukf.h"

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
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  //define augment spreading parameter
  lambda_aug = 3 - n_aug;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

    ///* Weights of sigma points
  //VectorXd weights_;
}

UKF::~UKF() {}

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

void AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

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
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_aug + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug + n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

void SigmaPointPrediction(double delta_t, MatrixXd Xsig_aug, MatrixXd* Xsig_out) {
  // predicted sigma points matrix
  Xsig_p = MatrixXd(n_x_, 2*n_aug_+1);
  //time diff in sec
  double delta_t = 0.1; 

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

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
    Xsig_p(0,i) = px_p;
    Xsig_p(1,i) = py_p;
    Xsig_p(2,i) = v_p;
    Xsig_p(3,i) = yaw_p;
    Xsig_p(4,i) = yawd_p;
  }

  //write result
  *Xsig_out = Xsig_p;
}

void PredictMeanAndCovariance(MatrixXd Xsig_pred_, VectorXd* x_pred, MatrixXd* P_pred) {
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_aug/(lambda_aug+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_aug);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_pred = x;
  *P_pred = P;

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

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(MatrixXd* Xsig_aug);
  SigmaPointPrediction(double delta_t, MatrixXd Xsig_aug, MatrixXd* Xsig_pred_);
  PredictMeanAndCovariance(MatrixXd Xsig_pred_, VectorXd* x_, MatrixXd* P_);

}


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
}
