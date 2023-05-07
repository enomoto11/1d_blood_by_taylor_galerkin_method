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

#include <cstdio>
#include <iomanip>

int Countnumofline(std::string file)
{
  std::ifstream ifs(file); // ファイルオープン
  int numofnode = 0;       // 行の数を数える定数
  std::string str;         // 文字列指定
  while (std::getline(ifs, str))
  {
    numofnode++;
  }
  ifs.close();
  return numofnode;
}

int main()
{
  int NODE_NUM = Countnumofline("input/node_d.dat");
  int ELEMENT_NUM = Countnumofline("input/element_d.dat");

  const double L = 1.0e-02;                                    // tubeの長さ[m]
  const int M = 1200;                                          // 時間ステップ数[-]
  const double dt = 1.0e-03;                                   // 時間刻み[s]
  const double DELTA_X = L / ELEMENT_NUM;                      // 要素の長さ[m]
  const double PI = M_PI;                                      // 円周率
  const double r0 = 1.0e-05;                                   // 初期状態のtubeの半径[m]
  const double h0 = 1.0e-06;                                   // 初期状態のtubeの厚さ[m]
  const double A0 = PI * pow(r0 - h0, 2.0e0);                  // 初期状態のtubeの流路面積（位置座標によらない）[m^2]
  const double K_R = 3.3e-06;                                  // 粘性抵抗K_R[m^2/s]
  const double rho = 1060000.0;                                // 密度[g/m^3]
  const double E = 1.0e08;                                     // ヤング率(0.1MPa)[g/m/s^2]
  const double betha = 4.0e0 / 3.0e0 * sqrt(PI) * h0 * E / A0; // ß
  // const double dbetha_dx = -4.0 / 3.0 * h0 * sqrt(PI) / PI / PI * E * 1.0e8; // dß/dx

  double v0 = 1.0e-2; // 所定位置における初期状態のtubeの流速[m/s]

  std::vector<std::vector<double>> node(NODE_NUM, std::vector<double>(3, 0));
  std::ifstream ifs("input/node_d.dat"); // ファイル入力
  for (int i = 0; i < NODE_NUM; i++)
  {
    std::vector<double> node_t(3, 0);
    for (int j = 0; j < 3; j++)
    {
      ifs >> node_t[j];
    }
    node[i] = node_t;
  }
  ifs.close();

  std::vector<std::vector<double>> element(ELEMENT_NUM, std::vector<double>(3, 0));
  std::fstream ifs1("input/element_d.dat"); // ファイル入力
  for (int i = 0; i < ELEMENT_NUM; i++)
  {
    std::vector<double> element_t(2, 0);
    for (int j = 0; j < 2; j++)
    {
      ifs1 >> element_t[j];
    }
    element[i] = element_t;
  }
  ifs1.close();

  using namespace Eigen;

  // [時刻][node番号]
  std::vector<std::vector<double>>
      flowQuantity(M + 1, std::vector<double>(NODE_NUM, 0.0));
  std::vector<std::vector<double>> area(M + 1, std::vector<double>(NODE_NUM, 0.0));
  std::vector<std::vector<double>> velocity(M + 1, std::vector<double>(NODE_NUM, 0.0));

  // init Q and A (which we want to calculate)
  for (int i = 0; i < NODE_NUM; i++)
  {
    area[0][i] = A0;
    // 初期状態のnodeに与える速度はノード点の位置によって変える
    if (i >= 0 && i <= NODE_NUM / 2)
    {
      velocity[0][i] = v0;
    }
    flowQuantity[0][i] = area[0][i] * velocity[0][i];
  }

  // 時間ステップごとに計算
  for (int i = 0; i < M; i++)
  {
    Eigen::MatrixXd A_area = MatrixXd::Zero(NODE_NUM, NODE_NUM);
    Eigen::VectorXd b_area = VectorXd::Zero(NODE_NUM);
    Eigen::MatrixXd A_flowQuantity = MatrixXd::Zero(NODE_NUM, NODE_NUM);
    Eigen::VectorXd b_flowQuantity = VectorXd::Zero(NODE_NUM);

    for (int j = 0; j < ELEMENT_NUM; j++)
    {
      int ele0 = element[j][0];
      int ele1 = element[j][1];

      A_area(ele0, ele0) = A_area(ele0, ele0) + 2.0e0 * DELTA_X / (dt * 6.0e0);
      A_area(ele0, ele1) = A_area(ele0, ele1) + 1.0e0 * DELTA_X / (dt * 6.0e0);
      A_area(ele1, ele0) = A_area(ele1, ele0) + 1.0e0 * DELTA_X / (dt * 6.0e0);
      A_area(ele1, ele1) = A_area(ele1, ele1) + 2.0e0 * DELTA_X / (dt * 6.0e0);

      A_flowQuantity(ele0, ele0) = A_flowQuantity(ele0, ele0) + 2.0e0 * DELTA_X / (dt * 6.0e0);
      A_flowQuantity(ele0, ele1) = A_flowQuantity(ele0, ele1) + 1.0e0 * DELTA_X / (dt * 6.0e0);
      A_flowQuantity(ele1, ele0) = A_flowQuantity(ele1, ele0) + 1.0e0 * DELTA_X / (dt * 6.0e0);
      A_flowQuantity(ele1, ele1) = A_flowQuantity(ele1, ele1) + 2.0e0 * DELTA_X / (dt * 6.0e0);
    }
    for (int j = 0; j < ELEMENT_NUM; j++)
    {
      int ele0 = element[j][0];
      int ele1 = element[j][1];
      double area1 = area[i][ele0], area2 = area[i][ele1];
      double flowQuantity1 = flowQuantity[i][ele0], flowQuantity2 = flowQuantity[i][ele1];

      // calculate first term in right side
      b_area(ele0) = b_area(ele0) + DELTA_X / 6.0e0 * (2.0e0 * area1 + 1.0e0 * area2);
      b_area(ele1) = b_area(ele1) + DELTA_X / 6.0e0 * (1.0e0 * area1 + 2.0e0 * area2);
      b_flowQuantity(ele0) = b_flowQuantity(ele0) + DELTA_X / 6.0e0 * (2.0e0 * flowQuantity1 + 1.0e0 * flowQuantity2);
      b_flowQuantity(ele1) = b_flowQuantity(ele1) + DELTA_X / 6.0e0 * (1.0e0 * flowQuantity1 + 2.0e0 * flowQuantity2);

      // calculate second term in right side
      b_area(ele0) = b_area(ele0) + dt * (-flowQuantity1 + flowQuantity2) / 2.0e0 - dt / 2.0e0 * K_R * (-(flowQuantity1 / area1) + (flowQuantity2 / area2) / 2.0e0);
      b_area(ele1) = b_area(ele1) + dt * (flowQuantity1 - flowQuantity2) / 2.0e0 - dt / 2.0e0 * K_R * ((flowQuantity1 / area1) - (flowQuantity2 / area2) / 2.0e0);
      b_flowQuantity(ele0) = b_flowQuantity(ele0) + dt * ((-(pow(flowQuantity1, 2.0e0) / area1) + pow(flowQuantity2, 2.0e0) / area2) / 2.0e0 + betha / rho / 3.0e0 / 2.0e0 * (-pow(area1, 1.5e0) + pow(area2, 1.5e0)) - dt * K_R / 2.0e0 * (-(pow(flowQuantity1 / area1, 2.0e0)) + pow(flowQuantity2 / area2, 2.0e0)));
      b_flowQuantity(ele1) = b_flowQuantity(ele1) + dt * ((pow(flowQuantity1, 2.0e0) / area1 - pow(flowQuantity2, 2.0e0) / area2) / 2.0e0 + betha / rho / 3.0e0 / 2.0e0 * (pow(area1, 1.5e0) - pow(area2, 1.5e0)) - dt * K_R / 2.0e0 * (pow(flowQuantity1 / area1, 2.0e0) - pow(flowQuantity2 / area2, 2.0e0)));

      // calculate third term in right side

      // calculate fourth term in right side

      // calculate fifth term in right side
      // b_area has no fifth term
      b_flowQuantity(ele0) = b_flowQuantity(ele0) - dt * (K_R * DELTA_X / 6.0e0 * (2.0e0 * (flowQuantity1 / area1) + 1.0e0 * (flowQuantity2 / area2)) + dt / 2.0e0 * pow(K_R, 2.0e0) * DELTA_X / 6.0e0 * (2.0e0 * (flowQuantity1 / pow(area1, 2.0e0)) + 1.0e0 * (flowQuantity2 / pow(area2, 2.0e0))));
      b_flowQuantity(ele1) = b_flowQuantity(ele1) - dt * (K_R * DELTA_X / 6.0e0 * (1.0e0 * (flowQuantity1 / area1) + 2.0e0 * (flowQuantity2 / area2)) + dt / 2.0e0 * pow(K_R, 2.0e0) * DELTA_X / 6.0e0 * (1.0e0 * (flowQuantity1 / pow(area1, 2.0e0)) + 2.0e0 * (flowQuantity2 / pow(area2, 2.0e0))));
    }

    // for (int j = 0; j < ELEMENT_NUM + 1; j++)
    // {
    //   A_area(0, j) = 0.0e0;
    //   A_flowQuantity(0, j) = 0.0e0;
    // }
    // A_area(0, 0) = 1.0e0;
    // A_flowQuantity(0, 0) = 1.0e0;
    // b_area(0) = 0.0e0;
    // b_flowQuantity(0) = 0.0e0;

    // 行列計算
    Eigen::VectorXd x_area = A_area.fullPivLu().solve(b_area);
    Eigen::VectorXd x_flowQuantity = A_flowQuantity.fullPivLu().solve(b_flowQuantity);

    for (int j = 0; j < x_area.size(); j++)
    {
      area[i + 1][j] = x_area(j);
    }
    for (int j = 0; j < x_flowQuantity.size(); j++)
    {
      flowQuantity[i + 1][j] = x_flowQuantity(j);
    }

    // TODO: debugを関数化したい
    // --------------------------------------------------------------------debug
    // b_areaの出力
    std::ostringstream oss_b_area;
    oss_b_area << "debug/b/area/debug" << i << ".csv";
    std::ofstream outputfile_b_area(oss_b_area.str());
    for (int j = 0; j < ELEMENT_NUM + 1; j++)
    {
      outputfile_b_area << b_area(j) << std::endl;
    }
    outputfile_b_area.close();

    // b_flowQuantityの出力
    std::ostringstream oss_b_flowQuantity;
    oss_b_flowQuantity << "debug/b/flowQuantity/debug" << i << ".csv";
    std::ofstream outputfile_b_flowQuantity(oss_b_flowQuantity.str());
    for (int j = 0; j < ELEMENT_NUM + 1; j++)
    {
      outputfile_b_flowQuantity << b_flowQuantity(j) << std::endl;
    }
    outputfile_b_flowQuantity.close();

    // x_areaの出力
    std::ostringstream oss_x_area;
    oss_x_area << "debug/x/area/debug" << i << ".csv";
    std::ofstream ofs_x_area(oss_x_area.str());
    for (int k = 0; k < x_area.size(); k++)
    {
      ofs_x_area << x_area(k) << std::endl;
    }
    ofs_x_area.close();

    // x_flowQuantityの出力
    std::ostringstream oss_x_flowQuantity;
    oss_x_flowQuantity << "debug/x/flowQuantity/debug" << i << ".csv";
    std::ofstream ofs_x_flowQuantity(oss_x_flowQuantity.str());
    for (int k = 0; k < x_flowQuantity.size(); k++)
    {
      ofs_x_flowQuantity << x_flowQuantity(k) << std::endl;
    }
    ofs_x_flowQuantity.close();

    // --------------------------------------------------------------------debug
  }

  int d = 100;
  // 見たい時刻t=0~200
  std::ofstream ofs("output/flowQuantity.dat");
  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < d + 1; j++)
    {
      ofs << flowQuantity[i * 200][j] << " ";
    }
    ofs << std::endl;
  }
  ofs.close();

  std::ofstream ofs2("output/area.dat");
  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < d + 1; j++)
    {
      ofs2 << area[i * 200][j] << " ";
    }
    ofs2 << std::endl;
  }

  return 0;
}