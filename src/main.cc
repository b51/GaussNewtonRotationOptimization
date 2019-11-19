/*************************************************************************
*
*              Author: b51
*                Mail: b51live@gmail.com
*            FileName: main.cc
*
*          Created On: Tue 19 Nov 2019 09:56:26 AM CST
*     Licensed under The MIT License [see LICENSE for details]
*
************************************************************************/

#include <iostream>
#include <fstream>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <Eigen/Core>
#include <vector>

#include "Undistortion.h"
#include "GaussNewton.h"

void ReadPointsFromFile(const std::string& file_name,
                        std::vector<Eigen::Vector3d>& Pws,
                        std::vector<Eigen::Vector2d>& ps) {
  Pws.clear();
  ps.clear();
  std::ifstream f(file_name);
  if (!f) {
    LOG(FATAL) << "The matchingPoints file does not exit!";
  } else {
    std::string str;
    while (getline(f, str)) {
      std::istringstream is(str);
      std::string s;
      std::vector<float> data;
      while (is >> s) {
        data.push_back(stof(s));
      }
      Pws.emplace_back(Eigen::Vector3d(data[0], data[1], data[2]));
      ps.emplace_back(Eigen::Vector2d(data[3], data[4]));
    }
  }
  f.close();
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  FLAGS_minloglevel = google::INFO;
  FLAGS_logtostderr = true;
  FLAGS_colorlogtostderr = true;

  std::vector<Eigen::Vector3d> Pws;
  std::vector<Eigen::Vector2d> ps;
  ReadPointsFromFile("../data/points.txt", Pws, ps);

  std::vector<Eigen::Vector2d> undistorted_ps;

  Eigen::Matrix3d K;
  Eigen::Matrix<double, 5, 1> dist_coeff;
  K << 579.1870, 0., 639.8485,
       0., 581.8611, 364.8964,
       0.,       0.,       1.;
  dist_coeff << -0.2866, 0.0954, 8.3960e-04, 0.0011, -0.0163;

  Undistortion undistortion(K, dist_coeff);
  undistortion.UndistortPoints(ps, undistorted_ps);

  Eigen::Vector3d t(0.0, 1.3299, -1.7923);
  GaussNewton gauss_newton(K);
  gauss_newton.SetTranslation(t);
  gauss_newton.Iteration(Pws, undistorted_ps, 100);
}
