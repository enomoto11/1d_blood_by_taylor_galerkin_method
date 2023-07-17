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

void FLOW1D::init(int toggle)
{
  if (toggle == 0)
  {
    for (int i = 0; i < NODE_NUM; i++)
    {
      area[i] = A0;
      velocity[i] = v0;
      flowQuantity[i] = area[i] * velocity[i];
    }
  }
  else if (toggle == 1)
  {
    for (int i = 0; i < NODE_NUM; i++)
    {
      double rate = 1.1e0;
      if (i == 0)
      {
        area[i] = A0;
        velocity[i] = 0e0;
        flowQuantity[i] = area[i] * velocity[i];
      }
      else if (i != 0 && i < NODE_NUM / 3)
      {
        area[i] = A0 / rate;
        velocity[i] = v0 * rate;
        flowQuantity[i] = area[i] * velocity[i];
      }
      else
      {
        area[i] = A0;
        velocity[i] = v0;
        flowQuantity[i] = area[i] * velocity[i];
      }
    }
  }
  else if (toggle == 2)
  {
    for (int i = 0; i < NODE_NUM; i++)
    {
      area[i] = A0;
      if (i == 0)
      {
        velocity[i + 1] = velocity[i];
        flowQuantity[i] = area[i] * velocity[i];
      }
      if (i == NODE_NUM - 2)
      {
        velocity[i] = velocity[i + 1];
      }
    }
  }
  else if (toggle == 3)
  {
    for (int i = 0; i < NODE_NUM; i++)
    {
      if (i < NODE_NUM / 2e0)
      {
        area[i] = A0;
        velocity[i] = v0;
        flowQuantity[i] = area[i] * velocity[i];
      }
      else
      {
        area[i] = A0;
        velocity[i] = v0 / 2e0;
        flowQuantity[i] = area[i] * velocity[i];
      }
    }
  }
  else if (toggle == 4)
  {
    for (int i = 0; i < NODE_NUM; i++)
    {
      area[i] = A0;
      velocity[i] = v0 * 10;
      flowQuantity[i] = area[i] * velocity[i];
    }
  }
  else if (toggle == 5)
  {
    area[0] = A0;
    velocity[0] = v0;
    flowQuantity[0] = area[0] * velocity[0];
    for (int i = 1; i < NODE_NUM; i++)
    {
      area[i] = A0;
      velocity[i] = v0 / 1e2;
      flowQuantity[i] = area[i] * velocity[i];
    }
  }
}
