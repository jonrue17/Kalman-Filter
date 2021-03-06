#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0; 
   
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     */

    // first measurement, initialise all the matrix 
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    ekf_.P_ = MatrixXd(4, 4); 
    ekf_.P_ << 0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0, 
               0, 0, 0, 0; 
    
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << 0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
        cout << "EKF: Radar Initialisation" << endl;
        ekf_.x_(0) = measurement_pack.raw_measurements_(0) / (sqrt(1.0 + pow(tan(measurement_pack.raw_measurements_(1)),2))); 
        ekf_.x_(1) = ekf_.x_(0) * tan(measurement_pack.raw_measurements_(1));
        float eps = 10 ^ -5; 
        if (abs(ekf_.x_(0)) > eps) { // avoid 0-division
            ekf_.x_(2) = measurement_pack.raw_measurements_(0) * measurement_pack.raw_measurements_(2) / ekf_.x_(0);
        }
        else if (abs(ekf_.x_(1)) > eps) { // avoid 0-division
           ekf_.x_(3) = measurement_pack.raw_measurements_(0) * measurement_pack.raw_measurements_(2) / ekf_.x_(1);
        }
        else {
            // keep 0m/s as initial speed 
        }
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state, use x and y position from LASER to initialise Px and Py State 
      // Whitout any speed measurement, set Vx and Vy to 0 m/s 
        ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0; 
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds. */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2) = dt; 
  ekf_.F_(1, 3) = dt;


   /* Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float noise_ax = process_noise_default; 
  float noise_ay = process_noise_default;

  ekf_.Q_(0, 0) = pow(dt, 4) / 4 * pow(noise_ax, 2); 
  ekf_.Q_(0, 2) = pow(dt, 3) / 2 * pow(noise_ax, 2);
  ekf_.Q_(1, 1) = pow(dt, 4) / 4 * pow(noise_ay, 2);
  ekf_.Q_(1, 3) = pow(dt, 3) / 2 * pow(noise_ay, 2);
  ekf_.Q_(2, 0) = pow(dt, 3) / 2 * pow(noise_ax, 2);
  ekf_.Q_(2, 2) = pow(dt, 2) * pow(noise_ax, 2);
  ekf_.Q_(3, 1) = pow(dt, 3) / 2 * pow(noise_ay, 2);
  ekf_.Q_(3, 3) = pow(dt, 2) * pow(noise_ay, 2);

  ekf_.P_ = ekf_.F_* ekf_.P_*(ekf_.F_.transpose())+ekf_.Q_;
  
  ekf_.Predict();
  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  Tools tools; 

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      cout << "EKF: Update Radar " << endl;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = MatrixXd(3, 4);
      ekf_.H_ = Hj_; 
      ekf_.R_ = MatrixXd(3, 3);
      ekf_.R_ = R_radar_; 
      ekf_.UpdateEKF(measurement_pack.raw_measurements_,ekf_.x_);

  } else {
    // Laser updates
      cout << "EKF: Update Laser " << endl;
      ekf_.H_ = MatrixXd(2, 4);
      ekf_.H_ = H_laser_; 
      ekf_.R_ = MatrixXd(2, 2);
      ekf_.R_ = R_laser_; 
      ekf_.Update(measurement_pack.raw_measurements_); 
  }
  cout << "EKF: Update Done " << endl;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
