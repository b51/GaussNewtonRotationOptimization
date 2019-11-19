/*************************************************************************
*
*              Author: b51
*                Mail: b51live@gmail.com
*            FileName: Undistortion.cc
*
*          Created On: Tue 19 Nov 2019 03:19:47 PM CST
*     Licensed under The MIT License [see LICENSE for details]
*
************************************************************************/

#include "Undistortion.h"

constexpr int N_MAX_ITERATION = 10;
constexpr double EPSILON = 1.;

Undistortion::Undistortion(const Eigen::Matrix3d& K,
                           const Eigen::Matrix<double, 5, 1> dist_coeff)
    : K_(K), dist_coeff_(dist_coeff) {}

Undistortion::~Undistortion() {}

void Undistortion::DistortPoint(const Eigen::Vector2d& p,
                                const Eigen::Matrix<float, 12, 1>& k,
                                Eigen::Vector2d& distorted_p) {
  double r2 = p[0] * p[0] + p[1] * p[1];
  double r4 = r2 * r2;
  double r6 = r4 * r2;
  double a1 = 2 * p[0] * p[1];
  double a2 = r2 + 2 * p[0] * p[0];
  double a3 = r2 + 2 * p[1] * p[1];
  double cdist = 1 + k[0] * r2 + k[1] * r4 + k[4] * r6;
  double icdist2 = 1. / (1 + k[5] * r2 + k[6] * r4 + k[7] * r6);
  distorted_p[0] =
      p[0] * cdist * icdist2 + k[2] * a1 + k[3] * a2 + k[8] * r2 + k[9] * r4;
  distorted_p[1] =
      p[1] * cdist * icdist2 + k[2] * a3 + k[3] * a1 + k[10] * r2 + k[11] * r4;
}

void Undistortion::UndistortPoint(const Eigen::Vector2d& p,
                                  Eigen::Vector2d& undistorted_p) {
  int dist_param_size = dist_coeff_.rows();
  Eigen::Matrix<float, 12, 1> k = Eigen::Matrix<float, 12, 1>::Zero();
  for (int i = 0; i < std::min(dist_param_size, 12); i++) k[i] = dist_coeff_[i];

  double fu = K_(0, 0);
  double fv = K_(1, 1);
  double cu = K_(0, 2);
  double cv = K_(1, 2);

  double x = (p[0] - cu) / fu;
  double y = (p[1] - cv) / fv;
  double x0 = x;
  double y0 = y;

  for (int i = 0; i < N_MAX_ITERATION; i++) {
    double r2 = x * x + y * y;
    double icdist = (1 + ((k[7] * r2 + k[6]) * r2 + k[5]) * r2) /
                    (1 + ((k[4] * r2 + k[1]) * r2 + k[0]) * r2);
    double deltaX =
        2 * k[2] * x * y + k[3] * (r2 + 2 * x * x) + k[8] * r2 + k[9] * r2 * r2;
    double deltaY = k[2] * (r2 + 2 * y * y) + 2 * k[3] * x * y + k[10] * r2 +
                    k[11] * r2 * r2;
    x = (x0 - deltaX) * icdist;
    y = (y0 - deltaY) * icdist;

    // distort point and compute error
    Eigen::Vector2d distorted_p;
    DistortPoint(Eigen::Vector2d(x, y), k, distorted_p);

    double x_proj = distorted_p[0] * fu + cu;
    double y_proj = distorted_p[1] * fv + cv;
    double error =
        std::sqrt(std::pow(x_proj - p[0], 2) + std::pow(y_proj - p[1], 2));

    if (error < EPSILON) break;
  }
  undistorted_p[0] = (int)(x * fu + cu + 0.5);
  undistorted_p[1] = (int)(y * fv + cv + 0.5);
}

void Undistortion::UndistortPoints(
    const std::vector<Eigen::Vector2d>& points,
    std::vector<Eigen::Vector2d>& undistorted_points) {
  size_t n = points.size();
  undistorted_points.clear();
  undistorted_points.reserve(n);
  for (const auto& p : points) {
    Eigen::Vector2d undistorted_p;
    UndistortPoint(p, undistorted_p);
    undistorted_points.emplace_back(undistorted_p);
  }
}
