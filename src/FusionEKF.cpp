#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

namespace {

const float NOISE_X = 9;
const float NOISE_Y = 9;

}
/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  ekf_.Rl_ = MatrixXd(2, 2);
  ekf_.Rr_ = MatrixXd(3, 3);
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  ekf_.Rl_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  ekf_.Rr_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.H_ << 1, 0, 0, 0,
            0, 1, 0, 0;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
  
  
  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
  
  
  }
/** * Destructor.  */ 
FusionEKF::~FusionEKF() {} 

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

/*****************************************************************************
*  Initialization
****************************************************************************/
if (!is_initialized_) {
/**
TODO:
* Initialize the state ekf_.x_ with the first measurement.
* Create the covariance matrix.
* Remember: you'll need to convert radar from polar to cartesian coordinates.
*/
// first measurement
cout << "EKF: " << endl;
VectorXd ekf_.x_(4); 

if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
/**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float fi = measurement_pack.raw_measurements_[1];
      float rdot = measurement_pack.raw_measurements_[2];

      ekf_.x_ << rho * cos(fi), rho * sin(fi), rdot * cos(fi), rdot * sin(fi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          1, 0, 1, 0,
          0, 1, 0, 1;
  ekf_.Q_(0, 0) = (noise_ax * pow(dt,4) / 4);
  ekf_.Q_(0, 2) = (noise_ax * pow(dt,3) /2);
  ekf_.Q_(1, 1) = (noise_ay * pow(dt,4) / 4);
  ekf_.Q_(1, 3) = (noise_ay * pow(dt,3) / 2);
  ekf_.Q_(2, 0) = (noise_ax * pow(dt,3) /2);
  ekf_.Q_(2, 2) = (noise_ax * pow(dt,2));
  ekf_.Q_(3, 1) = (noise_ay * pow(dt,3) / 2);
  ekf_.Q_(3, 3) = (noise_ay * pow(dt,2));

  ekf_.Predict();

/*****************************************************************************
*  Update
****************************************************************************/

/**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

}
