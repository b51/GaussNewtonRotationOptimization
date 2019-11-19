/*************************************************************************
 *
 *              Author: b51
 *                Mail: b51live@gmail.com
 *            FileName: GaussNewton.h
 *
 *          Created On: Tue 19 Nov 2019 10:24:19 AM CST
 *     Licensed under The MIT License [see LICENSE for details]
 *
 ************************************************************************/

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

class GaussNewton {
 public:
  GaussNewton(const Eigen::Matrix3d& K);
  ~GaussNewton();

  void Iteration(const std::vector<Eigen::Vector3d>& Pws,
                 const std::vector<Eigen::Vector2d>& ps, int n_iterations);

  void SetTranslation(const Eigen::Vector3d& t) { t_ = t; }

 private:
  Eigen::Matrix3d GetSkewMatrix(const Eigen::Vector3d& phi) {
    Eigen::Matrix3d skew;
    double a1 = phi[0];
    double a2 = phi[1];
    double a3 = phi[2];
    skew << 0., -a3, a2,
            a3, 0., -a1,
            -a2, a1, 0.;
    return skew;
  }

  Eigen::Matrix3d LieAlgebraToRotation(const Eigen::Vector3d& phi);

  Eigen::Matrix<double, 2, 3> JacobianCalculation(const Eigen::Matrix3d& R,
                                                  const Eigen::Vector3d& P);

  void IterateOnce(const std::vector<Eigen::Vector3d>& Pws,
                   const std::vector<Eigen::Vector2d>& ps);

  Eigen::Matrix<double, 2, 3> jacobian_;
  Eigen::Matrix3d K_;
  Eigen::Matrix3d R_;
  Eigen::Vector3d t_;
  int n_max_iteration_;
};
