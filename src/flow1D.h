#ifndef _FLOW1D_H_
#define _FLOW1D_H_

//##################################################################################
//
// 1D flow simulator
//
// Copyright (c) 2023 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   flow1D.h
 * @brief  Header
 * @author T. Otani, T. Enomoto
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdio>
#include <iomanip>

#include "shapefunction.h"
#include "gauss.h"
#include "right-part.h"
#include "exception-args.h"
#include "exception-manager.h"

#include "Eigen/Core"
#include "Eigen/LU"

class FLOW1D {
 public:
  

  int ELEMENT_NUM = 100;
  int NODE_NUM = ELEMENT_NUM + 1;

  std::vector<std::vector<double>> element = std::vector(ELEMENT_NUM, std::vector<double>(3, 0));
  std::vector<double> x;

  const double L = 1.0e-1;                                     // tubeの長さ[m]
  const int M = 1200;                                          // 時間ステップ数[-]
  const double dt = 2.0e-08;                                   // 時間刻み[s]
  const double DELTA_X = L / ELEMENT_NUM;                      // 要素の長さ[m]
  const double PI = M_PI;                                      // 円周率
  const double r0 = 1.0e-02;                                   // 初期状態のtubeの半径[m]
  const double h0 = r0 / 1e1;                                  // 初期状態のtubeの厚さ[m]
  const double A0 = PI * pow(r0 - h0, 2.0e0);                  // 初期状態のtubeの流路面積（位置座標によらない）[m^2]
  const double K_R = 3.3e0;                                    // 粘性抵抗K_R[m^2/s]
  const double rho = 1.06e6;                                   // 密度[g/m^3]
  const double E = 1.0e08;                                     // ヤング率(0.1MPa)[g/m/s^2]
  const double betha = 4.0e0 / 3.0e0 * sqrt(PI) * h0 * E / A0; // ß
  double v0 = 1.0e0;                                           // 所定位置における初期状態のtubeの流速[m/s]

  std::vector<double> area = std::vector(NODE_NUM, 0e0);
  std::vector<double> velocity = std::vector(NODE_NUM, 0e0);
  std::vector<double> flowQuantity = std::vector(NODE_NUM, 0e0);
  std::vector<double> pressure = std::vector(NODE_NUM, 0e0);

  void input();
  void init(int toggle);

  void output_init();
  void output(const int iter);
  void exec(const int iter);

  void exportVTP(const int iter);

 private:
  void compute_LHS(Eigen::MatrixXd &A);
  void compute_RHS(Eigen::VectorXd &b_area,Eigen::VectorXd &b_flowQuantity,const int iter);

};


#endif