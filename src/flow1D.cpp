#include "flow1D.h"

using namespace std;
using namespace Eigen;

void FLOW1D::exec(const int iter)
{
  MatrixXd A_area = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  MatrixXd A_flowQuantity = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  VectorXd b_area = VectorXd::Zero(NODE_NUM);
  VectorXd b_flowQuantity = VectorXd::Zero(NODE_NUM);

  compute_LHS(A_area);

  A_flowQuantity = A_area;

  compute_RHS(b_area,b_flowQuantity,iter);

  //boundary condition
  area[NODE_NUM - 1] = A0;

  Eigen::VectorXd x_area = A_area.fullPivLu().solve(b_area);

  //boundary condition
  for(int j=0;j<NODE_NUM;j++) A_flowQuantity(0,j) = 0e0;
  A_flowQuantity(0,0) = 1e0;
  flowQuantity[0] = A0 * v0 * (1e0 + 1e-1 * sin(2e0 * M_PI * 5e0 * iter / M));

  Eigen::VectorXd x_flowQuantity = A_flowQuantity.fullPivLu().solve(b_flowQuantity);

  for (int j = 0; j < NODE_NUM; j++)
  {
    area[j] = x_area(j);
    flowQuantity[j] = x_flowQuantity(j);
  }
}

void FLOW1D::compute_LHS(Eigen::MatrixXd &A)
{
  for (int j = 0; j < ELEMENT_NUM; j++)
  {
    int ele0 = element[j][0];
    int ele1 = element[j][1];

    Gauss g(1);
    vector<double> N(2);
    vector<double> dNdr(2);

    for (int k = 0; k < 2; k++)
    {
      ShapeFunction1D::P2_N(N, g.point[k]);
      ShapeFunction1D::P2_dNdr(dNdr, g.point[k]);

      double dxdr = dNdr.at(0) * x.at(ele0) + dNdr.at(1) * x.at(ele1);

      A(ele0, ele0) += N.at(0) * N.at(0) * g.weight[k] * dxdr;
      A(ele0, ele1) += N.at(0) * N.at(1) * g.weight[k] * dxdr;
      A(ele1, ele0) += N.at(1) * N.at(0) * g.weight[k] * dxdr;
      A(ele1, ele1) += N.at(1) * N.at(1) * g.weight[k] * dxdr;
    }
  }
}

void FLOW1D::compute_RHS(Eigen::VectorXd &b_area,Eigen::VectorXd &b_flowQuantity,const int iter)
{
  for (int j = 0; j < ELEMENT_NUM; j++)
  {
    int ele0 = element[j][0];
    int ele1 = element[j][1];
    double flow0 = flowQuantity[ele0], flow1 = flowQuantity[ele1];
    double area0 = area[ele0], area1 = area[ele1];

    Gauss g(1);
    vector<double> N(2);
    vector<double> dNdr(2);
    vector<double> dNdx(2);
    vector<double> dNinvdr(2);
    vector<double> dNinvdx(2);

    for (int k = 0; k < 2; k++)
    {
      ShapeFunction1D::P2_N(N, g.point[k]);
      ShapeFunction1D::P2_dNdr(dNdr, g.point[k]);
      ShapeFunction1D::P2_dNinvdr(dNinvdr, g.point[k]);

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

      //first term
      b_area(ele0) += N.at(0) * A * g.weight[k] * dxdr;
      b_area(ele1) += N.at(1) * A * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += N.at(0) * Q * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += N.at(1) * Q * g.weight[k] * dxdr;

      //second term
      b_area(ele0) += dNdx.at(0) * dt * (Q + dt / 2e0 * (-K_R * Q / A)) * g.weight[k] * dxdr;
      b_area(ele1) += dNdx.at(1) * dt * (Q + dt / 2e0 * (-K_R * Q / A)) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += dNdx.at(0) * dt * (Q * Q / A + betha / 3e0 / rho * pow(A, 1.5e0) + dt * Q / A * (-K_R * Q / A)) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += dNdx.at(1) * dt * (Q * Q / A + betha / 3e0 / rho * pow(A, 1.5e0) + dt * Q / A * (-K_R * Q / A)) * g.weight[k] * dxdr;

      //third term
      b_flowQuantity(ele0) += -N.at(0) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx - (K_R / A) * dQQAdx - (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += -N.at(1) * dt * dt / 2.0e0 * (K_R * Q / (A * A) * dQdx - (K_R / A) * dQQAdx - (betha * K_R / (2e0 * rho)) * pow(A, -5e-1) * dAdx) * g.weight[k] * dxdr;

      //fourth term
      b_area(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, 5e-1) * dAdx) * g.weight[k] * dxdr;
      b_area(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (dQQAdx + (betha / (2e0 * rho)) * pow(A, 5e-1) * dAdx) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += -dNdx.at(0) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += -dNdx.at(1) * dt * dt / 2.0e0 * (-(Q * Q / A) * dQdx + (betha / 2e0 / rho) * pow(A, 5e-1) * dQdx + 2e0 * Q / A * dQQAdx + betha / rho * pow(A, -5e-1) * Q * dAdx) * g.weight[k] * dxdr;

      //fifth term
      b_flowQuantity(ele0) += N.at(0) * dt * (-K_R * Q / A + dt / 2e0 * K_R * K_R * Q / A / A) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += N.at(1) * dt * (-K_R * Q / A + dt / 2e0 * K_R * K_R * Q / A / A) * g.weight[k] * dxdr;
    }
  }

}

void FLOW1D::output_init()
{
  int iteratedTime;

  fstream ifsIteratedTime("output/dat/iteratedTime.dat");
  ifsIteratedTime.close();

  ofstream ofs("output/dat/flowQuantity.dat");
  ofs.close();

  ofstream ofs2("output/dat/area.dat");
  ofs2.close();

  ofstream ofs3("output/dat/velocity.dat");
  ofs3.close();

  ofstream ofs4("output/dat/pressure1.dat");
  ofstream ofs5("output/dat/pressure2.dat");
  ofstream ofs6("output/dat/pressure3.dat");
  ofs4.close();
  ofs5.close();
  ofs6.close();
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
  ofs << endl;
  ofs.close();

  ofstream ofs2("output/dat/area.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs2 << area[j] << " ";
  }
  ofs2 << endl;
  ofs2.close();

  ofstream ofs3("output/dat/velocity.dat",ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs3 << flowQuantity[j] / area[j] << " ";
  }
  ofs3 << endl;
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

  ofs4 << endl;
  ofs5 << endl;
  ofs6 << endl;
  ofs4.close();
  ofs5.close();
  ofs6.close();
}


void FLOW1D::exportVTP(const int iter)
{

  string file = "output/vtp/data_" + to_string(iter) + ".vtp";

  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<PolyData>\n");

  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",NODE_NUM,1);
  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e 0e0 0e0\n",x[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Lines>\n");

  fprintf(fp, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%d ",i);
  }
  fprintf(fp,"\n");
  fprintf(fp,"</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  fprintf(fp,"%d\n",NODE_NUM);
  fprintf(fp,"</DataArray>\n");

  fprintf(fp, "</Lines>\n");

  fprintf(fp,"<CellData>\n");
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"LV0\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(int) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"parent\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(int) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"doughtor0\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(int) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"doughtor1\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(int) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(int) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(double) * numOfVessel;
  // fprintf(fp,"<DataArray type=\"Float64\" Name=\"flowMass\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // numOfData += sizeof(int) + sizeof(double) * numOfVessel;
  // if(vesselType==VesselType::ARTERY){
  // fprintf(fp,"<DataArray type=\"Int32\" Name=\"WillisID\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",numOfData);
  // }
  fprintf(fp,"</CellData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</PolyData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}