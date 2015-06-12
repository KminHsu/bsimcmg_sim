#include <iostream>
#include <math.h>
#include "dlfcn.h"
#include "mex.h"

using namespace std;

void BSIMCMG(char* psInstName, double Vd, double Vg, double Vs, double Ve);

typedef void (*BSIMCMGIF)(char*, double, double, double, double); 

void mexFunction(
		 int          nlhs,
		 mxArray      *[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
  double      *vin1, *vin2;

  /* Check for proper number of arguments */

  //if (nrhs != 2) {
  //  mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
  //          "MEXCPP requires two input arguments.");
  //} else if (nlhs >= 1) {
  //  mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
  //          "MEXCPP requires no output argument.");
  //}
 
  //void* lib_handle = dlopen("/mnt/hgfs/kminhsu/Desktop/Xyce-Build/Xyce-6.2-build/src/DeviceModelPKG/test/libbsimcmg_module.so", RTLD_LAZY); 
  //cout << "error:" << dlerror() << endl;
  //cout << lib_handle << endl; 
  //BSIMCMGIF BSIMCMG;
  //BSIMCMGCreateInst("M1", "nmos1", "L=3e-8") = (BSIMCMGIF) dlsym(lib_handle, "BSIMCMG");
  //if( NULL !=  BSIMCMG )
  //  BSIMCMG(1.0, 0.3, 0.0, 0.0 );
  //else  
  //  cout << "Can not get BSIMCMG function from plugin : " << dlerror() << endl;
 
  int ret = -1;
  size_t nBufLen = 0;
  char* psInstName = NULL;
  nBufLen = mxGetN(prhs[0])*sizeof(mxChar)+1;
  psInstName = (char*)mxMalloc(nBufLen);
  ret = mxGetString(prhs[0], psInstName, nBufLen);
 
  double Vd = mxGetScalar(prhs[1]);
  double Vg = mxGetScalar(prhs[2]);
  double Vs = mxGetScalar(prhs[3]);
  double Ve = mxGetScalar(prhs[4]);

  BSIMCMG(psInstName, Vd, Vg, Vs, Ve);

  //cout << psInstName << endl
  //     << Vd << endl
  //     << Vg << endl
  //     << Vs << endl
  //     << Ve << endl;
  return;
}
