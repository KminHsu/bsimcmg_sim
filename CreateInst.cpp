#include <iostream>
#include <math.h>
#include "dlfcn.h"
#include "mex.h"

#include <N_DEV_ADMSbsimcmg.h>

using namespace std;

Xyce::Device::ADMSbsimcmg::Instance* CreateInst(char* psInstName, char* pModelName, char* psInstParams);

typedef Xyce::Device::ADMSbsimcmg::Instance* (*CreateInstance)(char*, char*, char*);

void mexFunction(
		 int          nlhs,
		 mxArray      *[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{

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
  //CreateInstIF CreateInst;
  //CreateInst = (CreateInstIF) dlsym(lib_handle, "CreateInst");
  //if( NULL !=  CreateInst )
  //  CreateInst(1.0, 0.3, 0.0, 0.0 );
  //else  
  //  cout << "Can not get CreateInst function from plugin : " << dlerror() << endl;
 
  char* psModelName = NULL; 
  char* psInstParams = NULL; 

  int ret = -1;
  size_t nBufLen = 0;
  char* psInstName = NULL;
  nBufLen = mxGetN(prhs[0])*sizeof(mxChar)+1;
  psInstName = (char*)mxMalloc(nBufLen);
  ret = mxGetString(prhs[0], psInstName, nBufLen);

  nBufLen = mxGetN(prhs[1])*sizeof(mxChar)+1;
  psModelName = (char*)mxMalloc(nBufLen);
  ret = mxGetString(prhs[1], psModelName, nBufLen);
  
  nBufLen = mxGetN(prhs[2])*sizeof(mxChar)+1;
  psInstParams = (char*)mxMalloc(nBufLen);
  ret = mxGetString(prhs[2], psInstParams, nBufLen);

  //cout << "here:" << psInstName << " " << psModelName << " " << psInstParams << endl;
  Xyce::Device::ADMSbsimcmg::Instance* pInst = CreateInst(psInstName, psModelName, psInstParams);

  mxFree(psInstName);
  mxFree(psModelName);
  mxFree(psInstParams);
  return;
}
