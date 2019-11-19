/*************************************************************************
 *
 *              Author: b51
 *                Mail: b51live@gmail.com
 *            FileName: Undistortion.h
 *
 *          Created On: Tue 19 Nov 2019 03:16:43 PM CST
 *     Licensed under The MIT License [see LICENSE for details]
 *
 ************************************************************************/

#include <Eigen/Core>
#include <iostream>
#include <vector>

class Undistortion {
 public:
  Undistortion(const Eigen::Matrix3d& K,
               const Eigen::Matrix<double, 5, 1> dist_coeff);

  ~Undistortion();

  void UndistortPoints(const std::vector<Eigen::Vector2d>& points,
                       std::vector<Eigen::Vector2d>& undistorted_points);

 private:
  Eigen::Matrix3d K_;
  Eigen::Matrix<double, 5, 1> dist_coeff_;

  void DistortPoint(const Eigen::Vector2d& p,
                    const Eigen::Matrix<float, 12, 1>& k,
                    Eigen::Vector2d& distorted_p);

  void UndistortPoint(const Eigen::Vector2d& p, Eigen::Vector2d& undistorted_p);
};
