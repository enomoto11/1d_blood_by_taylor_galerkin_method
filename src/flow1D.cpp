#include "flow1D.h"

using namespace std;
using namespace Eigen;

void FLOW1D::exec(const int iter)
{
  MatrixXd A_area = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  VectorXd b_area = VectorXd::Zero(NODE_NUM);
  MatrixXd A_flowQuantity = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  VectorXd b_flowQuantity = VectorXd::Zero(NODE_NUM);

  compute_LHS(A_area,A_flowQuantity);

  // 境界条件
  area[NODE_NUM - 1] = A0;
  flowQuantity[0] = A0 * v0 * (1e0 + 1e-1 * sin(2e0 * M_PI * 5e0 * iter / M));

  compute_RHS(b_area,b_flowQuantity,iter);

  Eigen::VectorXd x_area = A_area.fullPivLu().solve(b_area);
  Eigen::VectorXd x_flowQuantity = A_flowQuantity.fullPivLu().solve(b_flowQuantity);

  for (int j = 0; j < NODE_NUM; j++)
  {
    area[j] = x_area(j);
    flowQuantity[j] = x_flowQuantity(j);
  }
}

void FLOW1D::compute_LHS(Eigen::MatrixXd &A_area,Eigen::MatrixXd &A_flowQuantity)
{
  ShapeFunction1D shape;
  for (int j = 0; j < ELEMENT_NUM; j++)
  {
    int ele0 = element[j][0];
    int ele1 = element[j][1];

    for (int k = 0; k < 2; k++)
    {
      Gauss g(1);
      vector<double> N(2, 0e0);
      vector<double> dNdr(2, 0e0);

      shape.P2_N(N, g.point[k]);
      shape.P2_dNdr(dNdr, g.point[k]);

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
}

void FLOW1D::compute_RHS(Eigen::VectorXd &b_area,Eigen::VectorXd &b_flowQuantity,const int iter)
{
  RightPart rightPart = RightPart::newRightPart();
  vector<int> indices = {}; // 右辺の1~5項の内、計算をskipするものを指定する, (example) {2,4}
  RightPart::disable(rightPart, indices);

  ShapeFunction1D shape;
  for (int j = 0; j < ELEMENT_NUM; j++)
  {
    int ele0 = element[j][0];
    int ele1 = element[j][1];
    double flow0 = flowQuantity[ele0], flow1 = flowQuantity[ele1];
    double area0 = area[ele0], area1 = area[ele1];

    for (int k = 0; k < 2; k++)
    {
      Gauss g(1);
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

      double Q = N.at(0) * flow0 + N.at(1) * flow1;
      double A = N.at(0) * area0 + N.at(1) * area1;
      double dQdx = dNdx.at(0) * flow0 + dNdx.at(1) * flow1;
      double dAdx = dNdx.at(0) * area0 + dNdx.at(1) * area1;
      double dAinvdx = dNinvdx.at(0) / area0 + dNinvdx.at(1) / area1;
      double dQQAdx = dAinvdx * Q * Q + (1 / A) * dQdx * Q + (1 / A) * Q * dQdx;

      if (rightPart.firstTerm.shouldCalculate)
      {
        b_area(ele0) += N.at(0) * A * g.weight[k] * dxdr;
        b_area(ele1) += N.at(1) * A * g.weight[k] * dxdr;
        b_flowQuantity(ele0) += N.at(0) * Q * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += N.at(1) * Q * g.weight[k] * dxdr;

        if (b_area(ele0) < 0e0 || b_area(ele1) < 0e0 || b_flowQuantity(ele0) < 0e0 || b_flowQuantity(ele1) < 0e0 || Q < 0e0 || A < 0e0)
        {
          ExceptionArgs args = ExceptionArgs::newCheckArgs(b_area(ele0), b_area(ele1), b_flowQuantity(ele0), b_flowQuantity(ele1), A, Q, iter, j, rightPart.firstTerm.index);
          ExceptionManager::check(args);
        }
      }
      if (rightPart.secondTerm.shouldCalculate)
      {
        b_area(ele0) += dNdx.at(0) * dt * (Q + dt / 2e0 * (-K_R * Q / A)) * g.weight[k] * dxdr;
        b_area(ele1) += dNdx.at(1) * dt * (Q + dt / 2e0 * (-K_R * Q / A)) * g.weight[k] * dxdr;
        b_flowQuantity(ele0) += dNdx.at(0) * dt * (Q * Q / A + betha / 3e0 / rho * pow(A, 1.5e0) + dt * Q / A * (-K_R * Q / A)) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += dNdx.at(1) * dt * (Q * Q / A + betha / 3e0 / rho * pow(A, 1.5e0) + dt * Q / A * (-K_R * Q / A)) * g.weight[k] * dxdr;

        if (b_area(ele0) < 0e0 || b_area(ele1) < 0e0 || b_flowQuantity(ele0) < 0e0 || b_flowQuantity(ele1) < 0e0 || Q < 0e0 || A < 0e0)
        {
          ExceptionArgs args = ExceptionArgs::newCheckArgs(b_area(ele0), b_area(ele1), b_flowQuantity(ele0), b_flowQuantity(ele1), A, Q, iter, j, rightPart.secondTerm.index);
          ExceptionManager::check(args);
        }
      }
      if (rightPart.thirdTerm.shouldCalculate)
      { // b_area has no third term
        b_flowQuantity(ele0) += -N.at(0) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx - (K_R / A) * dQQAdx - (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += -N.at(1) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx - (K_R / A) * dQQAdx - (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;

        if (b_area(ele0) < 0e0 || b_area(ele1) < 0e0 || b_flowQuantity(ele0) < 0e0 || b_flowQuantity(ele1) < 0e0 || Q < 0e0 || A < 0e0)
        {
          ExceptionArgs args = ExceptionArgs::newCheckArgs(b_area(ele0), b_area(ele1), b_flowQuantity(ele0), b_flowQuantity(ele1), A, Q, iter, j, rightPart.thirdTerm.index);
          ExceptionManager::check(args);
        }
      }
      if (rightPart.fourthTerm.shouldCalculate)
      {
        b_area(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, 5e-1) * dAdx) * g.weight[k] * dxdr;
        b_area(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, 5e-1) * dAdx) * g.weight[k] * dxdr;
        b_flowQuantity(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;

        if (b_area(ele0) < 0e0 || b_area(ele1) < 0e0 || b_flowQuantity(ele0) < 0e0 || b_flowQuantity(ele1) < 0e0 || Q < 0e0 || A < 0e0)
        {
          ExceptionArgs args = ExceptionArgs::newCheckArgs(b_area(ele0), b_area(ele1), b_flowQuantity(ele0), b_flowQuantity(ele1), A, Q, iter, j, rightPart.fourthTerm.index);
          ExceptionManager::check(args);
        }
      }
      if (rightPart.fifthTerm.shouldCalculate)
      { // b_area has no fifth term
        b_flowQuantity(ele0) += N.at(0) * dt * (-K_R * Q / A + dt / 2e0 * K_R * K_R * Q / A / A) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += N.at(1) * dt * (-K_R * Q / A + dt / 2e0 * K_R * K_R * Q / A / A) * g.weight[k] * dxdr;

        if (b_area(ele0) < 0e0 || b_area(ele1) < 0e0 || b_flowQuantity(ele0) < 0e0 || b_flowQuantity(ele1) < 0e0 || Q < 0e0 || A < 0e0)
        {
          ExceptionArgs args = ExceptionArgs::newCheckArgs(b_area(ele0), b_area(ele1), b_flowQuantity(ele0), b_flowQuantity(ele1), A, Q, iter, j, rightPart.fifthTerm.index);
          ExceptionManager::check(args);
        }
      }
    }
  }

}

void FLOW1D::output(const int iter)
{
  int iteratedTime;

  fstream ifsIteratedTime("output/dat/iteratedTime.dat",ios::app);
  ifsIteratedTime >> iteratedTime;
  ifsIteratedTime.close();

  ofstream ofs("output/dat/flowQuantity.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs << flowQuantity[j] << " ";
  }
  ofs.close();

  ofstream ofs2("output/dat/area.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs2 << area[j] << " ";
  }
  ofs2.close();

  ofstream ofs3("output/dat/velocity.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs3 << flowQuantity[j] / area[j] << " ";
  }
  ofs3.close();

  ofstream ofs4("output/dat/pressure1.dat",ios::app);
  ofstream ofs5("output/dat/pressure2.dat",ios::app);
  ofstream ofs6("output/dat/pressure3.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    if (j == NODE_NUM / 4)
    {
      ofs4 << betha * (sqrt(area[j]) - sqrt(A0)) << " ";
    }
    if (j == NODE_NUM / 2)
    {
      ofs5 << betha * (sqrt(area[j]) - sqrt(A0)) << " ";
    }
    if (j == NODE_NUM * 3 / 4)
    {
      ofs6 << betha * (sqrt(area[j]) - sqrt(A0)) << " ";
    }
  }
  ofs4.close();
  ofs5.close();
  ofs6.close();
}