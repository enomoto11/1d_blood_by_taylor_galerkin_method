#include "flow1D.h"

using namespace std;

void FLOW1D::input()
{
  fstream ifs1("../../input/element_d.dat");
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
  ifstream ifs2("../../input/node_d.dat");
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
}

void FLOW1D::init()
{
  const double area0 = getA0(0e0);
  const double velocity0 = initInflowVelocity(0);

  area[0] = area0;
  velocity[0] = velocity0;
  flowQuantity[0] = area[0] * velocity[0];
  for (int i = 1; i < NODE_NUM; i++)
  // for (int i = 0; i < NODE_NUM; i++)
  {
    area[i] = getA0(i / NODE_NUM);
    velocity[i] = velocity0 / 5e0;
    flowQuantity[i] = area[i] * velocity[i];
  }
}
