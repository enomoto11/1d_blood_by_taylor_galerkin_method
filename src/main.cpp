#include "flow1D.h"

#ifdef _GLOG
#include "glog/logging.h"
#endif

using namespace std;

int main(int argc, char *argv[])
{
#ifdef _GLOG
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif

  FLOW1D flow1D;

  flow1D.input();

  flow1D.init(5);
  flow1D.output_init();

  int output_iter = 10;

  // boundary condition
  flow1D.flowQuantity[0] = flow1D.v0 * flow1D.A0;

  double time = 0e0;

  flow1D.A_area = Eigen::MatrixXd::Zero(flow1D.NODE_NUM, flow1D.NODE_NUM);
  flow1D.A_flowQuantity = Eigen::MatrixXd::Zero(flow1D.NODE_NUM, flow1D.NODE_NUM);

  flow1D.compute_LHS(flow1D.A_area);
  flow1D.compute_LHS(flow1D.A_flowQuantity);

  for (int i = 0; i < flow1D.NODE_NUM; i++)
    flow1D.A_area(0, i) = 0e0;
  flow1D.A_area(0, 0) = 1e0;

  for (int i = 0; i < flow1D.NODE_NUM; i++)
    flow1D.A_flowQuantity(0, i) = 0e0;
  flow1D.A_flowQuantity(0, 0) = 1e0;

  for (int iter = 0; iter <= flow1D.iterMax; iter++)
  {
    time += flow1D.dt;

    flow1D.exec(iter);
    flow1D.output(iter);
    if (iter % output_iter == 0)
    {
      flow1D.exportVTP(iter / output_iter);
      flow1D.exportVTPWith3D(iter / output_iter);
    }
  }
}