#include "flow1D.h"

#ifdef _GLOG
#include "glog/logging.h"
#endif

using namespace std;

int main(int argc,char *argv[])
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

  //boundary condition
  double v0 = 1e-2;
  flow1D.flowQuantity[0] = v0 * flow1D.A0;

  double time = 0e0;

  for (int iter = 0; iter <= 3000; iter++)
  {
    time += flow1D.dt;

    flow1D.exec(iter);
    flow1D.output(iter);
    if(iter%output_iter==0){
      flow1D.exportVTP(iter/output_iter);
    }

  //boundary 
  flow1D.flowQuantity[0] = flow1D.A0 * v0;
  
  }
}