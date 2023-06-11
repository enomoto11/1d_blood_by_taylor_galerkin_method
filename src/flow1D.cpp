#include "flow1D.h"

using namespace std;
using namespace Eigen;

void FLOW1D::exec(const int iter)
{
  MatrixXd A_area = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  MatrixXd A_flowQuantity = MatrixXd::Zero(NODE_NUM, NODE_NUM);
  b_area = VectorXd::Zero(NODE_NUM);
  b_flowQuantity = VectorXd::Zero(NODE_NUM);

  compute_LHS(A_area);
  A_flowQuantity = A_area;

  compute_RHS(b_area, b_flowQuantity, iter);

  Eigen::VectorXd x_area = A_area.fullPivLu().solve(b_area);
  Eigen::VectorXd x_flowQuantity = A_flowQuantity.fullPivLu().solve(b_flowQuantity);

  for (int j = 0; j < NODE_NUM; j++)
  {
    area[j] = x_area(j);
    flowQuantity[j] = x_flowQuantity(j);
    pressure[j] = beta * (sqrt(area[j])-sqrt(A0));
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

void FLOW1D::compute_RHS(Eigen::VectorXd &b_area, Eigen::VectorXd &b_flowQuantity, const int iter)
{
  for (int j = 0; j < ELEMENT_NUM; j++)
  {
    int ele0 = element[j][0];
    int ele1 = element[j][1];
    double flow0 = flowQuantity[ele0], flow1 = flowQuantity[ele1];
    double area0 = area[ele0], area1 = area[ele1];

    Gauss g(1);
    vector<double> N(2),dNdr(2),dNdx(2);

    for (int k = 0; k < 2; k++)
    {
      ShapeFunction1D::P2_N(N, g.point[k]);
      ShapeFunction1D::P2_dNdr(dNdr, g.point[k]);

      double dxdr = dNdr[0] * x.at(ele0) + dNdr[1] * x.at(ele1);
      double drdx = 1e0 / dxdr;

      dNdx[0] = dNdr[0] * drdx;
      dNdx[1] = dNdr[1] * drdx;

      double Flow = N.at(0) * flow0 + N.at(1) * flow1;
      double Area = N.at(0) * area0 + N.at(1) * area1;

      double dFlow_dx = dNdx[0] * flow0 + dNdx[1] * flow1;
      double dAdx = dNdx[0] * area0 + dNdx[1] * area1;

      Vector2d G,B,dGdx;
      G(0) = Flow;
      G(1) = Flow*Flow/Area + beta/(3e0*rho)*pow(Area,1.5e0);

      dGdx(0) = dFlow_dx;
      dGdx(1) = 2e0*Flow/Area*dFlow_dx - Flow*Flow/(Area*Area)*dAdx + 5e-1*beta/rho*sqrt(Area)*dAdx;

      B(0) = 0e0;
      B(1) = -K_R*Flow/Area;

      Matrix2d dGdQ,dBdQ;
      dGdQ(0,0) = 0e0;
      dGdQ(0,1) = 1e0;
      dGdQ(1,0) = -Flow*Flow/(Area*Area)+5e-1*beta/rho*sqrt(Area);
      dGdQ(1,1) = 2e0*Flow/Area;

      dBdQ(0,0) = 0e0;
      dBdQ(0,1) = 0e0;
      dBdQ(1,0) = K_R*Flow/(Area*Area);
      dBdQ(1,1) =-K_R/Area;

      Vector2d G_LW = G + 5e-1*dt * dGdQ * B;
      Vector2d B_LW = B + 5e-1*dt * dBdQ * B;

      //first term
      b_area(ele0)         += N.at(0) * Area * g.weight[k] * dxdr;
      b_area(ele1)         += N.at(1) * Area * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += N.at(0) * Flow * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += N.at(1) * Flow * g.weight[k] * dxdr;

      //second term
      b_area(ele0)         += dNdx.at(0) * dt * G_LW(0) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += dNdx.at(0) * dt * G_LW(1) * g.weight[k] * dxdr;
      if(j!=ELEMENT_NUM-1){
        b_area(ele1)         += dNdx.at(1) * dt * G_LW(0) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += dNdx.at(1) * dt * G_LW(1) * g.weight[k] * dxdr;
      }

      //third term
      Vector2d tmp = dBdQ * dGdx;
      b_area(ele0)         += -N[0]*dt*dt/2e0 * tmp(0) * g.weight[k] * dxdr;
      b_area(ele1)         += -N[1]*dt*dt/2e0 * tmp(0) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += -N[0]*dt*dt/2e0 * tmp(1) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += -N[1]*dt*dt/2e0 * tmp(1) * g.weight[k] * dxdr;

      //fourth term
      tmp = dGdQ * dGdx;
      b_area(ele0)         += -dNdx.at(0) * dt * dt / 2e0 * tmp(0) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += -dNdx.at(0) * dt * dt / 2e0 * tmp(1) * g.weight[k] * dxdr;
      if(j!=ELEMENT_NUM-1){
        b_area(ele1)         += -dNdx.at(1) * dt * dt / 2e0 * tmp(0) * g.weight[k] * dxdr;
        b_flowQuantity(ele1) += -dNdx.at(1) * dt * dt / 2e0 * tmp(1) * g.weight[k] * dxdr;
      }

      //fifth term
      b_area(ele0)         += N.at(0) * dt * B_LW(0) * g.weight[k] * dxdr;
      b_area(ele1)         += N.at(1) * dt * B_LW(0) * g.weight[k] * dxdr;
      b_flowQuantity(ele0) += N.at(0) * dt * B_LW(1) * g.weight[k] * dxdr;
      b_flowQuantity(ele1) += N.at(1) * dt * B_LW(1) * g.weight[k] * dxdr;
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

  fstream ifsIteratedTime("output/dat/iteratedTime.dat", ios::app);
  ifsIteratedTime >> iteratedTime;
  ifsIteratedTime.close();

  ofstream ofs("output/dat/flowQuantity.dat", ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs << flowQuantity[j] << " ";
  }
  ofs << endl;
  ofs.close();

  ofstream ofs2("output/dat/area.dat", ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs2 << area[j] << " ";
  }
  ofs2 << endl;
  ofs2.close();

  ofstream ofs3("output/dat/velocity.dat", ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    ofs3 << flowQuantity[j] / area[j] << " ";
  }
  ofs3 << endl;
  ofs3.close();

  ofstream ofs4("output/dat/pressure1.dat", ios::app);
  ofstream ofs5("output/dat/pressure2.dat", ios::app);
  ofstream ofs6("output/dat/pressure3.dat", ios::app);
  for (int j = 0; j < NODE_NUM; j++)
  {
    if (j == NODE_NUM / 4)
    {
      ofs4 << beta * (sqrt(area[j]) - sqrt(A0)) << " ";
    }
    if (j == NODE_NUM / 2)
    {
      ofs5 << beta * (sqrt(area[j]) - sqrt(A0)) << " ";
    }
    if (j == NODE_NUM * 3 / 4)
    {
      ofs6 << beta * (sqrt(area[j]) - sqrt(A0)) << " ";
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

  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"Area\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",area[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",sqrt(area[i]/PI));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"FlowRate\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",flowQuantity[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",pressure[i]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float32\" Name=\"RHS_area\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",b_area(i));
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float32\" Name=\"RHS_Q\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<NODE_NUM;i++){
    fprintf(fp,"%e\n",b_flowQuantity[i]);
  }
  fprintf(fp,"</DataArray>\n");

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
  fprintf(fp,"</PointData>\n");

  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</PolyData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}