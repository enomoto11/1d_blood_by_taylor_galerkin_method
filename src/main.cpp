#include "flow1D.h"

using namespace std;

int main()
{
  FLOW1D flow1D;

  flow1D.input();

  flow1D.init(5);

  for (int iter = 0; iter < flow1D.M; iter++)
  {
    flow1D.exec(iter);
    flow1D.output(iter);
  }
}