#include <iostream>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Hj_in,
                        MatrixXd &H_in, MatrixXd &Rl_in, MatrixXd &Rr_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  Rr_ = Rr_in;
  Rl_ = Rl_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;    
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + Rl_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    x_ = x_ + K * y;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

    Tools::CalculateJacobian(x_, Hj_);
    
    VectorXd tmp(3);
    float sqrSumRt = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
    
    if (sqrSumRt == 0)
    {
        std::cout << "px == py == 0\n";
        return;
    }

    tmp(0) = sqrSumRt;
    tmp(1) = fabs(x_(0)) > 1e-4 ? atan2(x_(1), x_(0)) : 0;
    tmp(2) = fabs(sqrSumRt) > 1e-4 ? (x_(2) * x_(0) + x_(1) * x_(3)) / sqrSumRt : 0;
    
    VectorXd y = z - tmp;
    if (y(1) > M_PI)
    {
      y(1) -= 2 * M_PI;
    }

    if (y(1) < -1 * M_PI)
    {
      y(1) += 2 * M_PI;
    }

    MatrixXd S = Hj_ * P_ * Hj_.transpose() + Rr_;
    MatrixXd K = P_ * Hj_.transpose() * S.inverse();

    x_ = x_ + K * y;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I - K * Hj_) * P_;
}
