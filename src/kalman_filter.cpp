#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * Predict the state with linear function
   */
    x_ = F_ * x_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state by using Kalman Filter equations for linear measurement function (LASER)
   */
    // Calculate the error 
    VectorXd y(2); 
    y = z - H_ * x_; 

    // Calculate the optimal kalman gain 
    MatrixXd S_ = H_ * P_ * (H_.transpose()) + R_;
    MatrixXd Kopt_=P_* (H_.transpose()) *(S_.inverse());
    
    // Update with the optimal kalman gain
    x_ = x_ + Kopt_ * y; 
    P_ = P_ - Kopt_ * H_ * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const VectorXd& x_state) {
  /**
   * update the state by using Extended Kalman Filter equations for non-linear measurement function (RADAR)
   */
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // Link the measurment with state estimation with non-linear function  
    VectorXd h_x(3); 
    // Check division by zero
    if (fabs(px * px + py * py) < 0.0001) {
        h_x << sqrt(px * px + py * py),
            std::atan2(py,px), 
            0;
    }
    else {
        h_x << sqrt(px * px + py * py),
            std::atan2(py, px),
            (px * vx + py * vy) / sqrt(px * px + py * py);
    }
    
    // Calculate the error 
    VectorXd y(3);
    y = z - h_x;

    // Making sure -pi<Phi<pi
    while (y(1) > 3.1415) {
        y(1) = y(1) - 2*3.1415; 
    }
    while (y(1) < -3.1415) {
        y(1) = y(1) + 2*3.1415;
    }
    std::cout << "Phi = " << y(1) << std::endl; 

    // Calculate the optimal kalman gain 
    MatrixXd S_ = H_ * P_ * (H_.transpose()) + R_;
    MatrixXd Kopt_ = P_ * (H_.transpose()) * (S_.inverse());

    // Update with the optimal kalman gain
    x_ = x_ + Kopt_ * y;
    P_ = P_ - Kopt_ * H_ * P_;
}
