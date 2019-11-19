/*************************************************************************
 *
 *              Author: b51
 *                Mail: b51live@gmail.com
 *            FileName: GaussNewton.cc
 *
 *          Created On: Tue 19 Nov 2019 12:21:36 PM CST
 *     Licensed under The MIT License [see LICENSE for details]
 *
 ************************************************************************/

#include "GaussNewton.h"

#include <glog/logging.h>
#include <gflags/gflags.h>

constexpr double EPSILON = 1e-5;

GaussNewton::GaussNewton(const Eigen::Matrix3d& K) : K_(K) {
  R_ = Eigen::Matrix3d::Identity();
}

GaussNewton::~GaussNewton() {}

Eigen::Matrix3d GaussNewton::LieAlgebraToRotation(const Eigen::Vector3d& phi) {
  double theta = phi.norm();
  Eigen::Vector3d a = phi / theta;
  Eigen::Matrix3d skew_a = GetSkewMatrix(a);

  Eigen::Matrix3d R = cos(theta) * Eigen::Matrix3d::Identity() +
                      (1 - cos(theta)) * a * a.transpose() +
                      sin(theta) * skew_a;
  return R;
}

/**
 *  Jacobian Calculation with Rotaiton
 *  Reference: https://www.cnblogs.com/gaoxiang12/p/5689927.html
 *
 *  error(phi + d_phi) = error(phi) + K * d_phi^ * exp(phi^) * P
 *  令 q = d_phi^ * exp(phi^) * P;
 *  令 u = K * q;
 *  对 u 在 0 处进行泰勒展开
 *  ==>
 *  u(d_phi) = u(0) + (du/dq) * (dq/d(d_phi)) * d_phi
 *
 *  du / dq = [du/dX du/dY du/dZ
 *             dv/dX dv/dY dv/dZ]
 *
 *  dq/d(d_phi) = d(d_phi^ * exp(phi^) * P)/d(d_phi)
 *              = -(R*P)^*d(d_phi) / d(d_phi)
 *              = -(R*P)^
 *
 *  (du/dq) * (dq/d(d_phi)) = [-fu*X*Y/Z_2     fu+fu*X_2/Z_2  -fu*Y/Z
 *                             -fv-fv*Y_2/Z_2  fv*X*Y/Z_2      fv*X/Z]
 *
 *  e(phi+d_phi) = e(phi) + Jacobian * d_phi
 *  Jacobian = (du/dq) * (dq/d(d_phi))
 */
Eigen::Matrix<double, 2, 3> GaussNewton::JacobianCalculation(
    const Eigen::Matrix3d& R, const Eigen::Vector3d& P) {
  double fu = K_(0, 0);
  double fv = K_(1, 1);
  Eigen::Matrix<double, 2, 3> jacobian;
  Eigen::Vector3d RP = R * P;
  double X = RP[0];
  double Y = RP[1];
  double Z = RP[2];
  double X_2 = X * X;
  double Y_2 = Y * Y;
  double Z_2 = Z * Z;

  double j00 = -fu * X * Y / Z_2;
  double j01 = fu + fu * X_2 / Z_2;
  double j02 = -fu * Y / Z;
  double j10 = -fv - fv * Y_2 / Z_2;
  double j11 = fv * X * Y / Z_2;
  double j12 = fv * X / Z;

  jacobian << j00, j01, j02,
              j10, j11, j12;
  return -1. * jacobian;
}

void GaussNewton::IterateOnce(const std::vector<Eigen::Vector3d>& Pws,
                              const std::vector<Eigen::Vector2d>& ps) {
  CHECK_EQ(Pws.size(), ps.size());

  size_t n = Pws.size();
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();

  for (size_t i = 0; i < n; i++) {
    Eigen::Vector3d Pw = Pws[i];
    Eigen::Vector2d p = ps[i];
    Eigen::Matrix<double, 2, 3> jacobian = JacobianCalculation(R_, Pw);
    A += (jacobian.transpose() * jacobian);

    Eigen::Vector3d p_homo = K_ * (R_ * Pw + t_);
    Eigen::Vector2d p_projection = (p_homo / p_homo[2]).block<2, 1>(0, 0);

    Eigen::Vector2d error = p - p_projection;
    b += (-1. * jacobian.transpose() * error);
  }
  /**
   *  This function solve Ax=b with LU decomposition
   */
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
  Eigen::Vector3d delta_phi = lu.solve(b);
  Eigen::Matrix3d delta_R = LieAlgebraToRotation(delta_phi);

  R_ = delta_R * R_;
}

Eigen::Matrix3d GaussNewton::Iteration(const std::vector<Eigen::Vector3d>& Pws,
                                       const std::vector<Eigen::Vector2d>& ps,
                                       int n_iterations) {
  double last_eTe = 0.;
  for (int i = 0; i < n_iterations; i++) {
    IterateOnce(Pws, ps);
    double eTe = 0.;
    for (size_t j = 0; j < Pws.size(); j++) {
      Eigen::Vector3d Pw = Pws[j];
      Eigen::Vector2d p = ps[j];
      Eigen::Vector3d p_homo = K_ * (R_ * Pw + t_);
      Eigen::Vector2d p_projection = (p_homo / p_homo[2]).block<2, 1>(0, 0);
      Eigen::Vector2d error = p - p_projection;
      eTe += error.transpose() * error;
      LOG(INFO) << "points[" << j << "] error: " << error.transpose();
    }
    if (eTe < EPSILON or std::abs(eTe - last_eTe) < EPSILON) {
      LOG(INFO) << "Optimization done with " << i << " times iteration";
      LOG(INFO) << "error: " << eTe;
      LOG(INFO) << "Optimized Rotation: \n" << R_;
      return R_;
    }
    LOG(INFO) << "iteration: " << i;
    LOG(INFO) << "R: \n" << R_;
    LOG(INFO) << "error: " << eTe;
    last_eTe = eTe;
  }
}
