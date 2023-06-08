
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

  int output_iter = 1;

  for (int iter = 0; iter < 1; iter++)
  {
    flow1D.exec(iter);
    flow1D.output(iter);
    if(iter%output_iter==0){
      flow1D.exportVTP(iter/output_iter);
    }
  }
}