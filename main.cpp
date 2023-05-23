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
#include <math.h>
#include <cmath>

#include "lib/eigen/Eigen/Core"
#include "lib/eigen/Eigen/LU"

#include "shapefunction.h"
#include "gauss.h"
#include "flag.h"
#include "params.h"

#include <cstdio>
#include <iomanip>

using namespace std;
using namespace Eigen;

int exec(LeftPartFlag flag)
{
  ShapeFunction1D shape;

  vector<vector<double>> element(ELEMENT_NUM, vector<double>(3, 0));
  fstream ifs1("input/element_d.dat"); // ファイル入力
  for (int i = 0; i < ELEMENT_NUM; i++)
  {
    vector<double> element_t(2, 0);
    for (int j = 0; j < 2; j++)
    {
      ifs1 >> element_t[j];
    }
    element[i] = element_t;
  }
  ifs1.close();

  string str;
  ifstream ifs2("input/node_d.dat");
  vector<double> x;
  while (getline(ifs2, str))
  {
    istringstream ss(str);
    string tmp;
    vector<double> tmp_x;
    for (int j = 0; j < 1; j++)
    {
      getline(ss, tmp, ' ');
      tmp_x.push_back(stod(tmp));
    }
    x.push_back(tmp_x[0]);
  }
  ifs2.close();

  // init Q and A (which we want to calculate)
  for (int i = 0; i < NODE_NUM; i++)
  {
    // 初期状態のnodeに与える速度はノード点の位置によって変える
    area[0][i] = A0;
    velocity[0][i] = v0;
    flowQuantity[0][i] = area[0][i] * velocity[0][i];
  }

  // 時間ステップごとに計算
  for (int i = 0; i < M; i++)
  {
    MatrixXd A_area = MatrixXd::Zero(NODE_NUM, NODE_NUM);
    VectorXd b_area = VectorXd::Zero(NODE_NUM);
    MatrixXd A_flowQuantity = MatrixXd::Zero(NODE_NUM, NODE_NUM);
    VectorXd b_flowQuantity = VectorXd::Zero(NODE_NUM);

    for (int j = 0; j < ELEMENT_NUM; j++)
    {
      Gauss g(1);
      int ele0 = element[j][0];
      int ele1 = element[j][1];

      for (int k = 0; k < 2; k++)
      {
        vector<double> N(2, 0e0);
        vector<double> dNdr(2, 0e0);
        vector<double> dNdx(2, 0e0);
        vector<double> dNinvdr(2, 0e0);
        vector<double> dNinvdx(2, 0e0);

        shape.P2_N(N, g.point[k]);
        shape.P2_dNdr(dNdr, g.point[k]);
        shape.P2_dNinvdr(dNinvdr, g.point[k]);

        double dxdr = dNdr.at(0) * x.at(ele0) + dNdr.at(1) * x.at(ele1);

        A_area(ele0, ele0) += N.at(0) * N.at(0) * g.weight[k] * dxdr;
        A_area(ele0, ele1) += N.at(0) * N.at(1) * g.weight[k] * dxdr;
        A_area(ele1, ele0) += N.at(1) * N.at(0) * g.weight[k] * dxdr;
        A_area(ele1, ele1) += N.at(1) * N.at(1) * g.weight[k] * dxdr;

        A_flowQuantity(ele0, ele0) += N.at(0) * N.at(0) * g.weight[k] * dxdr;
        A_flowQuantity(ele0, ele1) += N.at(0) * N.at(1) * g.weight[k] * dxdr;
        A_flowQuantity(ele1, ele0) += N.at(1) * N.at(0) * g.weight[k] * dxdr;
        A_flowQuantity(ele1, ele1) += N.at(1) * N.at(1) * g.weight[k] * dxdr;
      }
    }
    for (int j = 0; j < ELEMENT_NUM; j++)
    {
      Gauss g(1);
      // B.C. : 0番目のnode点では常に流路面積,速度一定
      area[i][0] = A0;
      flowQuantity[i][0] = A0 * v0;

      int ele0 = element[j][0];
      int ele1 = element[j][1];
      double area1 = area[i][ele0], area2 = area[i][ele1];
      double flowQuantity1 = flowQuantity[i][ele0], flowQuantity2 = flowQuantity[i][ele1];

      // calculate first term in right side
      if (flag.shouldCalculateFirstTerm)
      {
        for (int k = 0; k < 2; k++)
        {
          vector<double> N(2, 0e0);
          vector<double> dNdr(2, 0e0);
          vector<double> dNdx(2, 0e0);
          vector<double> dNinvdr(2, 0e0);
          vector<double> dNinvdx(2, 0e0);

          shape.P2_N(N, g.point[k]);
          shape.P2_dNdr(dNdr, g.point[k]);
          shape.P2_dNinvdr(dNinvdr, g.point[k]);

          double dxdr = dNdr.at(0) * x.at(ele0) + dNdr.at(1) * x.at(ele1);

          double Q = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];
          double A = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];

          b_area(ele0) += N.at(0) * A * g.weight[k] * dxdr;
          b_area(ele1) += N.at(1) * A * g.weight[k] * dxdr;
          b_flowQuantity(ele0) += N.at(0) * Q * g.weight[k] * dxdr;
          b_flowQuantity(ele1) += N.at(1) * Q * g.weight[k] * dxdr;
        }
      }

      // calculate second term in right side
      if (flag.shouldCalculateSecondTerm)
      {
        b_area(ele0) = b_area(ele0) + dt * (-flowQuantity1 + flowQuantity2) / 2.0e0 - dt / 2.0e0 * K_R * (-(flowQuantity1 / area1) + (flowQuantity2 / area2) / 2.0e0);
        b_area(ele1) = b_area(ele1) + dt * (flowQuantity1 - flowQuantity2) / 2.0e0 - dt / 2.0e0 * K_R * ((flowQuantity1 / area1) - (flowQuantity2 / area2) / 2.0e0);
        b_flowQuantity(ele0) = b_flowQuantity(ele0) + dt * ((-(pow(flowQuantity1, 2.0e0) / area1) + pow(flowQuantity2, 2.0e0) / area2) / 2.0e0 + betha / rho / 3.0e0 / 2.0e0 * (-pow(area1, 1.5e0) + pow(area2, 1.5e0)) - dt * K_R / 2.0e0 * (-(pow(flowQuantity1 / area1, 2.0e0)) + pow(flowQuantity2 / area2, 2.0e0)));
        b_flowQuantity(ele1) = b_flowQuantity(ele1) + dt * ((pow(flowQuantity1, 2.0e0) / area1 - pow(flowQuantity2, 2.0e0) / area2) / 2.0e0 + betha / rho / 3.0e0 / 2.0e0 * (pow(area1, 1.5e0) - pow(area2, 1.5e0)) - dt * K_R / 2.0e0 * (pow(flowQuantity1 / area1, 2.0e0) - pow(flowQuantity2 / area2, 2.0e0)));
      }

      // calculate third term in right side
      if (flag.shouldCalculateThirdTerm)
      {
        for (int k = 0; k < 2; k++)
        {
          vector<double> N(2, 0e0);
          vector<double> dNdr(2, 0e0);
          vector<double> dNdx(2, 0e0);
          vector<double> dNinvdr(2, 0e0);
          vector<double> dNinvdx(2, 0e0);

          shape.P2_N(N, g.point[k]);
          shape.P2_dNdr(dNdr, g.point[k]);
          shape.P2_dNinvdr(dNinvdr, g.point[k]);

          double dxdr = dNdr.at(0) * x.at(ele0) + dNdr.at(1) * x.at(ele1);
          double drdx = 1e0 / dxdr;

          dNdx.at(0) = dNdr.at(0) * drdx;
          dNdx.at(1) = dNdr.at(1) * drdx;

          dNinvdx.at(0) = dNinvdr.at(0) * drdx;
          dNinvdx.at(1) = dNinvdr.at(1) * drdx;

          double Q = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];
          double A = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];
          double dQdx = dNdx.at(0) * flowQuantity[i][ele0] + dNdx.at(1) * flowQuantity[i][ele1];
          double dAdx = dNdx.at(0) * area[i][ele0] + dNdx.at(1) * area[i][ele1];
          double dAinvdx = dNinvdx.at(0) / area[i][ele0] + dNinvdx.at(1) / area[i][ele1];
          double dQQAdx = dAinvdx * Q * Q + (1 / A) * dQdx * Q + (1 / A) * Q * dQdx;

          // b_area has no third term
          b_flowQuantity(ele0) += -N.at(0) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx + (K_R / A) * dQQAdx + (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
          b_flowQuantity(ele1) += -N.at(1) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx + (K_R / A) * dQQAdx + (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
        }
      }

      // calculate fourth term in right side
      if (flag.shouldCalculateFourthTerm)
      {
        for (int k = 0; k < 2; k++)
        {
          vector<double> N(2, 0e0);
          vector<double> dNdr(2, 0e0);
          vector<double> dNdx(2, 0e0);
          vector<double> dNinvdr(2, 0e0);
          vector<double> dNinvdx(2, 0e0);

          shape.P2_N(N, g.point[k]);
          shape.P2_dNdr(dNdr, g.point[k]);
          shape.P2_dNinvdr(dNinvdr, g.point[k]);

          double dxdr = dNdr.at(0) * x.at(ele0) + dNdr.at(1) * x.at(ele1);
          double drdx = 1e0 / dxdr;

          dNdx.at(0) = dNdr.at(0) * drdx;
          dNdx.at(1) = dNdr.at(1) * drdx;

          dNinvdx.at(0) = dNinvdr.at(0) * drdx;
          dNinvdx.at(1) = dNinvdr.at(1) * drdx;

          double Q = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];
          double A = N.at(0) * flowQuantity[i][ele0] + N.at(1) * flowQuantity[i][ele1];
          double dQdx = dNdx.at(0) * flowQuantity[i][ele0] + dNdx.at(1) * flowQuantity[i][ele1];
          double dAdx = dNdx.at(0) * area[i][ele0] + dNdx.at(1) * area[i][ele1];
          double dAinvdx = dNinvdx.at(0) / area[i][ele0] + dNinvdx.at(1) / area[i][ele1];
          double dQQAdx = dAinvdx * Q * Q + (1 / A) * dQdx * Q + (1 / A) * Q * dQdx;

          b_area(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
          b_area(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
          b_flowQuantity(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;
          b_flowQuantity(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;
        }
      }

      // calculate fifth term in right side
      if (flag.shouldCalculateFifthTerm)
      {
        // b_area has no fifth term
        b_flowQuantity(ele0) = b_flowQuantity(ele0) - dt * (K_R * DELTA_X / 6.0e0 * (2.0e0 * (flowQuantity1 / area1) + 1.0e0 * (flowQuantity2 / area2)) + dt / 2.0e0 * pow(K_R, 2.0e0) * DELTA_X / 6.0e0 * (2.0e0 * (flowQuantity1 / pow(area1, 2.0e0)) + 1.0e0 * (flowQuantity2 / pow(area2, 2.0e0))));
        b_flowQuantity(ele1) = b_flowQuantity(ele1) - dt * (K_R * DELTA_X / 6.0e0 * (1.0e0 * (flowQuantity1 / area1) + 2.0e0 * (flowQuantity2 / area2)) + dt / 2.0e0 * pow(K_R, 2.0e0) * DELTA_X / 6.0e0 * (1.0e0 * (flowQuantity1 / pow(area1, 2.0e0)) + 2.0e0 * (flowQuantity2 / pow(area2, 2.0e0))));
      }
    }

    // 行列計算
    Eigen::VectorXd x_area = A_area.fullPivLu().solve(b_area);
    Eigen::VectorXd x_flowQuantity = A_flowQuantity.fullPivLu().solve(b_flowQuantity);

    for (int j = 0; j < NODE_NUM; j++)
    {
      area[i + 1][j] = x_area(j);
      flowQuantity[i + 1][j] = x_flowQuantity(j);
    }

    // TODO: debugを関数化したい
    // --------------------------------------------------------------------debug
    if (i < 10 || i % 50 == 0)
    {
      // b_areaの出力
      ostringstream oss_b_area;
      oss_b_area << "debug/b/area/debug" << i << ".csv";
      ofstream outputfile_b_area(oss_b_area.str());
      for (int j = 0; j < ELEMENT_NUM + 1; j++)
      {
        outputfile_b_area << b_area(j) << endl;
      }
      outputfile_b_area.close();

      // b_flowQuantityの出力
      ostringstream oss_b_flowQuantity;
      oss_b_flowQuantity << "debug/b/flowQuantity/debug" << i << ".csv";
      ofstream outputfile_b_flowQuantity(oss_b_flowQuantity.str());
      for (int j = 0; j < ELEMENT_NUM + 1; j++)
      {
        outputfile_b_flowQuantity << b_flowQuantity(j) << endl;
      }
      outputfile_b_flowQuantity.close();

      // x_areaの出力
      ostringstream oss_x_area;
      oss_x_area << "debug/x/area/debug" << i << ".csv";
      ofstream ofs_x_area(oss_x_area.str());
      for (int k = 0; k < x_area.size(); k++)
      {
        ofs_x_area << x_area(k) << endl;
      }
      ofs_x_area.close();

      // x_flowQuantityの出力
      ostringstream oss_x_flowQuantity;
      oss_x_flowQuantity << "debug/x/flowQuantity/debug" << i << ".csv";
      ofstream ofs_x_flowQuantity(oss_x_flowQuantity.str());
      for (int k = 0; k < x_flowQuantity.size(); k++)
      {
        ofs_x_flowQuantity << x_flowQuantity(k) << endl;
      }
      ofs_x_flowQuantity.close();
    }
    // --------------------------------------------------------------------debug
  }

  int d = 100;
  ofstream ofs("output/flowQuantity.dat");
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < d + 1; j++)
    {
      ofs << flowQuantity[i * 100][j] << " ";
    }
    ofs << endl;
  }
  ofs.close();

  ofstream ofs2("output/area.dat");
  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < d + 1; j++)
    {
      ofs2 << area[i * 100][j] << " ";
    }
    ofs2 << endl;
  }

  return 0;
}

int main()
{
  LeftPartFlag lpf;
  lpf.shouldCalculateSecondTerm = false;
  lpf.shouldCalculateFourthTerm = false;
  lpf.shouldCalculateFifthTerm = false;

  exec(lpf);
}