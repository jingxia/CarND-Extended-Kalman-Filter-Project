#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
          || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
  
    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
                rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

void Tools::CalculateJacobian(const VectorXd& x_state, MatrixXd& Hj) {
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  if (px == 0 && py == 0)
  {
      std::cout << "The px and py are 0 \n";
  }
  else
  {
      float sqrSum =pow(px, 2) + pow(py, 2);
      Hj(0, 0) = px / sqrt(sqrSum);
      Hj(0, 1) = py / sqrt(sqrSum);
      Hj(0, 2) = 0;
      Hj(0, 3) = 0;
      Hj(1, 0) = -1 * py / sqrSum;
      Hj(1, 1) = px / sqrSum;
      Hj(1, 2) = 0;
      Hj(1, 3) = 0;
      Hj(2, 0) = py * (vx * py - vy * px) / (sqrt(sqrSum) * sqrSum);
      Hj(2, 1) = px * (vy * px - vx * py) / (sqrt(sqrSum) * sqrSum);
      Hj(2, 2) = px / sqrt(sqrSum);
      Hj(2, 3) = py / sqrt(sqrSum);
  }
  
}
